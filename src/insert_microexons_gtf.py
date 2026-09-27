"""Add ME_DB microexons to a GTF, for tools that need them as annotated exons.

The input GTF is not changed; a new GTF is written. For each microexon that is
not already an exon (same chromosome, strand, start and end), one host
transcript is chosen: a transcript of the same strand with an intron that
contains the microexon, leaving at least one intron base on each side. Hosts
are ranked MANE_Select, Ensembl_canonical, basic, then more exons, then
transcript ID, so the choice is deterministic. A copy of the host's transcript
and exon lines, with the microexon added, is appended as a new transcript
`<host>_ME_<start>_<end>` tagged `microexon_insert`. CDS, UTR and codon lines
are not copied, since the frame may change.

One host per microexon is enough for Whippet, which can join any annotated
donor to any annotated acceptor within a gene, and it keeps the number of
extra transcripts small.

ME_DB lines may be BED12 (the inner blocks are microexons), BED6
(0-based start) or a single `chrom_strand_start_end` ID (0-based start),
mixed in one file. Output coordinates are GTF (1-based, closed).

Also writes a report: one row per microexon with its status (present,
inserted, no_host) and the host transcript.
"""

import argparse
import gzip
import re
from collections import defaultdict


def _open(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def read_me_db(path):
    """Set of (chrom, strand, start, end), 1-based closed."""
    microexons = set()
    with _open(path) as stream:
        for line in stream:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 12:
                start = int(fields[1])
                sizes = [int(x) for x in fields[10].strip(",").split(",")]
                offsets = [int(x) for x in fields[11].strip(",").split(",")]
                for size, offset in list(zip(sizes, offsets))[1:-1]:
                    microexons.add((fields[0], fields[5], start + offset + 1, start + offset + size))
            elif len(fields) >= 6:
                microexons.add((fields[0], fields[5], int(fields[1]) + 1, int(fields[2])))
            else:
                chrom, strand, start, end = fields[0].strip().rsplit("_", 3)
                microexons.add((chrom, strand, int(start) + 1, int(end)))
    return microexons


def _attribute(attributes, key):
    match = re.search(r'(?:^|;\s*){} "([^"]*)"'.format(key), attributes)
    return match.group(1) if match else None


def _rank(transcript):
    tags = transcript["tags"]
    return (0 if "MANE_Select" in tags else 1, 0 if "Ensembl_canonical" in tags else 1,
            0 if "basic" in tags else 1, -len(transcript["exons"]), transcript["id"])


def insert(gtf_path, microexons):
    """Returns (gtf text, report text)."""
    lines, transcripts, exons_seen = [], {}, set()
    with _open(gtf_path) as stream:
        for line in stream:
            lines.append(line if line.endswith("\n") else line + "\n")
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("transcript", "exon"):
                continue
            transcript_id = _attribute(fields[8], "transcript_id")
            if transcript_id is None:
                continue
            record = transcripts.setdefault(transcript_id, {
                "id": transcript_id, "chrom": fields[0], "strand": fields[6],
                "transcript_line": None, "exon_lines": [], "exons": [], "tags": set()})
            record["tags"].update(re.findall(r'tag "([^"]+)"', fields[8]))
            if fields[2] == "transcript":
                record["transcript_line"] = fields
            else:
                start, end = int(fields[3]), int(fields[4])
                record["exons"].append((start, end))
                record["exon_lines"].append(fields)
                exons_seen.add((fields[0], fields[6], start, end))

    introns = defaultdict(list)   # (chrom, strand) -> [(start, end, transcript_id)]
    for record in transcripts.values():
        blocks = sorted(record["exons"])
        for (_, left_end), (right_start, _) in zip(blocks, blocks[1:]):
            if right_start - left_end > 1:
                introns[(record["chrom"], record["strand"])].append(
                    (left_end + 1, right_start - 1, record["id"]))

    report = ["chrom\tstrand\tstart\tend\tstatus\thost\n"]
    added = []
    for chrom, strand, start, end in sorted(microexons):
        if (chrom, strand, start, end) in exons_seen:
            report.append("{}\t{}\t{}\t{}\tpresent\tNA\n".format(chrom, strand, start, end))
            continue
        hosts = [transcripts[tid] for i_start, i_end, tid in introns[(chrom, strand)]
                 if i_start < start and end < i_end]
        if not hosts:
            report.append("{}\t{}\t{}\t{}\tno_host\tNA\n".format(chrom, strand, start, end))
            continue
        host = min(hosts, key=_rank)
        new_id = "{}_ME_{}_{}".format(host["id"], start, end)
        template = host["exon_lines"][0]

        def retag(fields, feature, first, last):
            attributes = fields[8].replace('transcript_id "{}"'.format(host["id"]),
                                           'transcript_id "{}"'.format(new_id))
            attributes = re.sub(r'\s*exon_number "?[^";]*"?;', "", attributes)
            attributes = re.sub(r'\s*exon_id "[^"]*";', "", attributes)
            attributes = attributes.rstrip() + ' tag "microexon_insert";'
            return "\t".join(fields[:2] + [feature, str(first), str(last)] + fields[5:8] + [attributes]) + "\n"

        if host["transcript_line"] is not None:
            line = host["transcript_line"]
            added.append(retag(line, "transcript", line[3], line[4]))
        blocks = sorted(host["exons"] + [(start, end)], reverse=(strand == "-"))
        for number, (first, last) in enumerate(blocks, start=1):
            added.append(retag(template, "exon", first, last).rstrip("\n")
                         + ' exon_number "{}";\n'.format(number))
        report.append("{}\t{}\t{}\t{}\tinserted\t{}\n".format(chrom, strand, start, end, host["id"]))
    return "".join(lines + added), "".join(report)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--me-db", required=True)
    parser.add_argument("--output", required=True, help="GTF (.gz for gzip)")
    parser.add_argument("--report", required=True)
    args = parser.parse_args(argv)
    text, report = insert(args.gtf, read_me_db(args.me_db))
    if args.output.endswith(".gz"):
        with open(args.output, "wb") as stream:
            stream.write(gzip.compress(text.encode(), mtime=0))
    else:
        with open(args.output, "w") as stream:
            stream.write(text)
    with open(args.report, "w") as stream:
        stream.write(report)


if __name__ == "__main__":
    main()

"""Stream complete FASTQ records and verify ordered mate identifiers."""

import argparse
import gzip
import re
import sys
from itertools import zip_longest
from pathlib import Path


def open_fastq(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def records(path, run_id):
    with open_fastq(path) as stream:
        number = 0
        while True:
            header = stream.readline()
            if not header:
                return
            number += 1
            sequence = stream.readline()
            separator = stream.readline()
            quality = stream.readline()
            if not sequence or not separator or not quality:
                raise ValueError("{} record {}: incomplete four-line FASTQ in {}".format(run_id, number, path))
            if not header.startswith("@") or not separator.startswith("+"):
                raise ValueError("{} record {}: invalid FASTQ header or separator in {}".format(run_id, number, path))
            if len(sequence.rstrip("\r\n")) != len(quality.rstrip("\r\n")):
                raise ValueError("{} record {}: sequence/quality length mismatch in {}".format(run_id, number, path))
            yield header.rstrip("\r\n")


def mate_id(header):
    token = header[1:].split()[0]
    return re.sub(r"/[12]$", "", token)


def mate_number(header):
    fields = header[1:].split()
    suffix = re.search(r"/([12])$", fields[0])
    if suffix:
        return suffix.group(1)
    if len(fields) > 1 and re.match(r"^[12]:", fields[1]):
        return fields[1][0]
    return None


def validate(read1, read2=None, run_id="unknown"):
    if read2 is None:
        for _ in records(read1, run_id):
            pass
        return
    for number, pair in enumerate(zip_longest(records(read1, run_id), records(read2, run_id)), 1):
        first, second = pair
        if first is None or second is None:
            raise ValueError("{} record {}: unequal mate counts".format(run_id, number))
        if mate_id(first) != mate_id(second):
            raise ValueError("{} record {}: mate ID mismatch: {} versus {}".format(
                run_id, number, first, second))
        if mate_number(first) not in (None, "1") or mate_number(second) not in (None, "2"):
            raise ValueError("{} record {}: mate designation mismatch".format(run_id, number))
        first_fields, second_fields = first[1:].split(), second[1:].split()
        if len(first_fields) > 1 and len(second_fields) > 1:
            f1, f2 = first_fields[1], second_fields[1]
            if re.match(r"^[12]:", f1) and re.match(r"^[12]:", f2) and not (f1.startswith("1:") and f2.startswith("2:")):
                raise ValueError("{} record {}: Illumina mate designation mismatch".format(run_id, number))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("R1")
    parser.add_argument("R2", nargs="?")
    parser.add_argument("--marker", required=True)
    parser.add_argument("--run-id", default="unknown")
    args = parser.parse_args(argv)
    marker = Path(args.marker)
    try:
        validate(args.R1, args.R2, args.run_id)
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text("validated\n")
    except (ValueError, OSError, EOFError, UnicodeError) as error:
        marker.unlink(missing_ok=True)
        print(str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

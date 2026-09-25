"""Keep one splice-junction tag per distinct annotated flank path.

Tag flanks are cut from spliced transcripts. When transcripts that share an
exon boundary continue differently within a flank's length (another upstream
or downstream exon, or a different exon end), each path gives a different
flank sequence. Keeping only one of them loses every read from the other
isoforms that reaches past the point where the paths diverge. These helpers
reduce the annotated paths to the distinct, maximal ones: identical sequences
are merged, and a sequence that is fully contained in a longer one at the same
junction is dropped, because any read that fits it also fits the longer one.
No unannotated combinations are created.
"""


def maximal_flanks(sequences, side):
    """Distinct flanks, dropping those contained in a longer one.

    side "up": flanks end at the junction, so a flank is redundant if another
    one ends with it. side "down": flanks start at the junction, so a flank is
    redundant if another one starts with it.
    """
    distinct = sorted(set(sequences), key=lambda s: (-len(s), s))
    kept = []
    for seq in distinct:
        if side == "up":
            contained = any(other.endswith(seq) for other in kept)
        else:
            contained = any(other.startswith(seq) for other in kept)
        if not contained:
            kept.append(seq)
    return kept


def maximal_tags(tags):
    """Distinct (up_len, sequence, down_len, ...) tags with the junction at up_len.

    A tag is dropped if another tag contains it with the junctions aligned.
    Extra fields after the third are carried along unchanged; the first tag of
    each distinct sequence is kept.
    """
    ordered = sorted(tags, key=lambda t: (-(t[0] + t[2]), t[1]))
    kept = []
    for tag in ordered:
        up, seq, down = tag[0], tag[1], tag[2]
        contained = any(
            other[0] >= up and other[2] >= down
            and other[1][other[0] - up: other[0] + down] == seq
            for other in kept
        )
        if not contained:
            kept.append(tag)
    return kept


def path_label(transcript, index):
    """Tag header transcript field: unchanged for the first path, #pN after."""
    return transcript if index == 0 else "{}#p{}".format(transcript, index)

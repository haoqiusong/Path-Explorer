import sys

input_path = sys.argv[1]
output_path = sys.argv[2]

with open(input_path, "r", encoding="utf-8", errors="replace") as fin, open(output_path, "w") as fout:
    header = ""
    seq = ""
    for line in fin:
        line = line.strip().replace("\r", "")
        if line.startswith(">"):
            if header and seq:
                fout.write(header + "\n")
                fout.write(seq + "\n")
            header = line
            seq = ""
        else:
            seq += line.upper().replace("U", "T")
    if header and seq:
        fout.write(header + "\n")
        fout.write(seq + "\n")
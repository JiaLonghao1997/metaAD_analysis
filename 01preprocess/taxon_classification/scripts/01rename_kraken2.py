import sys

infile = sys.argv[1]
outfile = sys.argv[2]
splitstr = sys.argv[3]
namespace = int(sys.argv[4])

with open(infile) as taxons, open(outfile, "w") as output:
    for taxon in taxons:
        if taxon.startswith("0"):
            names = taxon.strip().split("\t")
            outline = ''
            for name in names:
                if name == '0':
                    outline = "clade_name"
                else:
                    newname = name.split('-bracken')[0].split(splitstr)[namespace]
                    outline = outline + '\t' + newname

            output.write(outline + "\n")
        else:
            # continue
            output.write(taxon)

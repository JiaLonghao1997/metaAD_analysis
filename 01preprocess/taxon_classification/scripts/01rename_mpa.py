import sys

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile) as taxons, open(outfile, "w") as output:
    for taxon in taxons:
        if taxon.startswith("clade_name"):
            names = taxon.strip().split("\t")
            outline = ''
            for name in names:
                if name == 'clade_name':
                    outline = name
                else:
                    newname = name.split('_')[1]
                    outline = outline + '\t' + newname

            output.write(outline + "\n")
        else:
            # continue
            output.write(taxon)

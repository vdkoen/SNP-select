
"""
Scripts takes as input a raw unfilterd vcf file, and a file containing the groups
it will use the groups and select only variants which are the same within the groups.
This is necessary if you have multiple samples of the same species within the group.
For making a robust snp set you dont want snps to be selected which can differ the same species.
"""

def main():

    somedict = {}
    include = []
    exclude = ['E3','G12','H1','H2','H3','H4','H5','H6','H7','H8','H9','H10', "B4", "C5 ", "C9", "C12", "D7", "E12", "F2", "F5", "F8", "F11", "G4", "G11"]

    print(len(exclude))

    with open("rasnaam_lijst.txt", "r") as file:
        for line in file:

            try:
                splitted_line = (line.split("\t"))

                if splitted_line[2] in somedict:
                    # append the new number to the existing array at this slot
                    somedict[splitted_line[2]].append(splitted_line[10])
                else:
                    # create a new array in this slot
                    somedict[splitted_line[2]] = [splitted_line[10]]
            except IndexError:
                pass
    print(somedict)

    # pairs = open("pairs.txt", 'w')
    for x in somedict.values():
        for y in x:
            if y not in exclude:
                include.append(y)



        #
        #     pairs.write(str(y) + "\t")
        # pairs.write("\n")
    print(','.join((include)))
    # pairs.close()

    with open("pairs.txt", "r") as file2:
        for line in file2:
            y = line.split()
            for x in y:
                if x not in include:
                    print(x)

main()
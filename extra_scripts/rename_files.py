# -*- coding: utf-8 -*-
"""renames all the code names in the files to the real name.
file names should look like this: code1_somthing.faztq
new name will look like this: name1_something.fastq

input format:
name1   code1
name2   code2
name3   code3
"""
import os
def main():

    somedict = {}
    name_list = []

    with open("grassen_names.txt", "r") as file:
        for line in file:
            splitted_line = line.split()
            try:
                if splitted_line[0] in somedict:
                    # append the new number to the existing array at this slot
                    somedict[splitted_line[0]].append(splitted_line[1])
                else:
                    # create a new array in this slot
                    somedict[splitted_line[0]] = [splitted_line[1]]
            except IndexError:
                pass

    print(somedict)
    teller = 0
    check = 0
    for filename in os.listdir("/nfs/BigData01/Big_Data/Lolium/KG000150/KG000150/fastq"):
        if not filename.startswith("pool"):
            # print(filename)
            splitted_name = filename.split("_")

            # print(splitted_name[0][4:])

            # splitted_x = splitted_name.split(".")




            # if (splitted_name[0][:-1]) != "Ace":
            #
            #     print(splitted_name[0], end=",")
            #     check += 1
            # if check == 5:
            #     print()
            #     check = 0

            for name, value in somedict.items():
                # print(splitted_name, value)
                if splitted_name[0] in value:
                    name_list.append(name)
                    counter = name_list.count(name)
                    splitted_name[0] = name + ".P_Pool"

                    teller += 1

                    print(filename, '_'.join(splitted_name))
                    #filename
                    command = "mv /nfs/BigData01/Big_Data/Lolium/raw_pools/%s /nfs/BigData01/Big_Data/Lolium/raw_pools/%s" % (filename, '.'.join(splitted_name))
                    print(command)
                    # os.system('echo "%s" | sudo -S %s' % ("your_sudo_password", command))





    # for key in somedict.keys():
    #     for x in range(0, name_list.count(key)):
    #         print(key + str(x+1), end=",")
    #     print()

        # print(key, name_list.count(key))


    # print(name_list)
# /nfs/BigData01/Big_Data/L15-217_framboos/data/raw_sequences_names
# framboos_names_pairs.txt

if __name__ == '__main__':
    main()

# sudoPassword = 'mypass'
# command = 'mount -t vboxsf myfolder /home/myuser/myfolder'
# p = os.system('echo %s|sudo -S %s' % (sudoPassword, command))

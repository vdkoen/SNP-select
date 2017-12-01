import os
import gzip
input_dir = '/Users/thomasg/Downloads/360_VCFs_2.50'
output_dir = '/tmp/360_VCFs_2.50'

output_file = '/tmp/merged.vcf.gz'

for f in os.listdir(input_dir):
    if f.startswith('.'):
        continue
    if os.path.exists(os.path.join(output_dir,f)) or os.path.exists(os.path.join(output_dir,f.replace('.gz',''))):
        continue
    output_file = open(os.path.join(output_dir,f.replace('.gz','')),'w')
    with gzip.open(os.path.join(input_dir,f)) as file_handle:
        for n,line in enumerate(file_handle):
            if n and not n%100000:
                print('%s lines processed of %s' % (n,f))
            line = line.decode()
            if line[0] == "#":
                output_file.write(line)
            else:
                split_line = line.split("\t")
                split_line[8] = 'GT:DP:AD'
                value_dict = {}
                try:
                    for k,v in [l.split('=') for l in split_line[7].split(';')]:
                        value_dict[k] = v
                        if k == 'DP4':
                            k2 = 'AD'
                            v2 = sum([int(i) for n,i in enumerate(v.split(',')) if n >=2])
                            value_dict[k2] = v2
                            value_dict['DPHQ'] = sum([int(i) for n,i in enumerate(v.split(','))])
                    split_line[9] = split_line[9].rstrip('\n')
                    split_line[9] += ':%(DPHQ)s:%(AD)s' % value_dict
                    output_file.write('\t'.join(split_line)+'\n')
                except ValueError:
                    pass

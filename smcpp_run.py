import os
import sys
import time

start_time = time.time()
vcfin = ''
outpath = ''
threads_num = '1'
mutation_R = '2.5e-08'
generation_time = '2'

for i in range(len(sys.argv)):
    if sys.argv[i] == '-vcf':
        vcfin = sys.argv[i + 1]
    if sys.argv[i] == '-o':
        outpath = sys.argv[i + 1]
    if sys.argv[i] == '-t':
        threads_num = sys.argv[i + 1]
    if sys.argv[i] == '-m':
        mutation_R = sys.argv[i + 1]
    if sys.argv[i] == '-g':
        generation_time = sys.argv[i + 1]
    if sys.argv[i] == '-h' or len(sys.argv) <= 2:
        print('''SMC++ pipeline
        #### Please run "conda activate smcpp" at first! ####
        Usage: python smcpp_run.py -vcf sample.vcf -o /out/put/folder -t 20

        -vcf    the VCF file from gatk4_pipeline.py or other callSNP pipeline
        -o      output prefix
        -t      threads number

        optional parameters:
        -m      mutation rate (default: 2.5e-08)
        -g      generation time (default: 2)
                Reference:
                    human:  25
                    rat  :  1
                    fish :  1~2 (under normal condition)

        ''')
        exit(0)

def readVCFInfo(vcfpath_in):
    all_genome_length = 0
    seq_length_list = []
    sample_name = []
    with open(vcfpath_in) as f:
        for line in f:
            if '#' in line:
                if '##contig=<ID=' in line:
                    seq_length_list.append([line.split('=')[2].split(',')[0], int(line.split('=')[3][:-2])])
                    all_genome_length += int(line.split('=')[3][:-2])
                if '#CHROM' in line:
                    sample_name += line.split()[9:]
            else:
                break
    if len(seq_length_list) == 0 or len(sample_name) == 0:
        print('###ERROR!! VCF file not include chromosome length and sample info, please check again.')
        exit(0)
    
    seq_length_list = sorted(seq_length_list, reverse=True, key=lambda a:a[1])
    for i in range(len(seq_length_list)):
        # print('chr' + i, seq_length_list[i])
        if float(seq_length_list[i][1]) / all_genome_length <= 0.01:
            seq_length_list = seq_length_list[:i + 1]
            break

    return seq_length_list, ','.join(sample_name)


if __name__ == '__main__':
    if 'usage: smc++ [-h]' not in os.popen('smc++ -h').read():
        print('###ERROR!! Please run \"conda activate smcpp\" to change conda environment at first! ###')
    seq_name_list, sampleIDs = readVCFInfo(vcfin)
    os.system(f'mkdir -p {outpath}')
    print(f'{outpath} Start!')
    print('1. Compress & index VCF file')
    bgzip_outname = outpath + '/' + vcfin.split('/')[-1] + '.gz'
    bgzip_realname = vcfin.split('/')[-1] + '.gz'
    outname = outpath.split('/')[-1]
    os.system(f'bgzip -c {vcfin} > {bgzip_outname}')
    os.system(f'tabix {bgzip_outname}')
    os.chdir(outpath)
    log_out = '1> smc.log 2>&1'
    print('2. Change vcf to smc format')
    count_chr = 0
    all_chr_num = str(len(seq_name_list))
    for s in seq_name_list:
        count_chr += 1
        chr_name = str(count_chr)
        print(f'{chr_name}/{all_chr_num} \nRunning command: smc++ vcf2smc {bgzip_realname} chr_{s[0]}.smc.gz {s[0]} {outname}:{sampleIDs}')
        os.system(f'smc++ vcf2smc {bgzip_realname} chr_{s[0]}.smc.gz {s[0]} {outname}:{sampleIDs} {log_out}')
        print('\t', s[0], 'finished')
    print('3. Estimate Population Hisory')
    print(f'Running command: smc++ estimate --spline cubic --knots 15 --timepoints 1000 10000000 --cores {threads_num} -o {outname} {mutation_R} *.smc.gz')
    os.system(f'smc++ estimate --spline cubic --knots 15 --timepoints 1000 10000000 --cores {threads_num} -o {outname} {mutation_R} *.smc.gz {log_out}')
    print('\tEstimate finished')
    print(f'4. Plot result to {outname}.pdf')
    os.system(f'smc++ plot {outname}.pdf {outname}/model.final.json -g {generation_time} -c {log_out}')
    end_time = time.time()
    time_for_h = str(float(end_time - start_time) / 3600)[:4]
    print(f'All Processed ! Use {time_for_h} hours\n')

    
    

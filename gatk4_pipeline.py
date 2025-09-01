###################################################
#
#   gatk4_pipeline.py
#
#   This script is for automatic gatk4.1.8.1 single sample gvcf to joint call vcf pipeline
#
#   Written by Lu Liang
#
#   2025.6.17 at SHOU
#
###################################################

import os
import time
import sys
from multiprocessing import Pool
import copy
import json


referenceGenome = ''
fqFilePath = ''
outputPath = ''
t_num = '1'
old_db = ''

for i in range(len(sys.argv)):
    if sys.argv[i] == '-ref':
        referenceGenome = sys.argv[i + 1]
    if sys.argv[i] == '-fq':
        fqFilePath = sys.argv[i + 1]
    if sys.argv[i] == '-o':
        outputPath = sys.argv[i + 1]
    if sys.argv[i] == '-t':
        t_num = sys.argv[i + 1]
    if sys.argv[i] == '-db':
        old_db = sys.argv[i + 1]
    if sys.argv[i] == '-h' or sys.argv[i] == '--help' or len(sys.argv) <= 2:
        print('''
This script is for automatic gatk4.1.8.1 pipeline, single sample call then joint call.

Usage: 
conda activate gatk4
python3 gatk4_pipeline.py -ref genome.fa -fq /fq/file/path/ -o /output/path/ -t 10
        
Request Software:   bwa
                    samtools
                    gatk4
                    java8

    Necessary parameters:
        -ref            reference genome(genome.fa), use softmasked genome better.
        -fq             fq or fq.gz file folder, will use all .fq/.fastq/fq.gz/fastq.gz file in this folder
                        Make sure file name like xxx_R1.fq/xxx_R2.fq or xxx_R1.fq.gz/xxx_R2.fq.gz
        -o              all files output path
                        # You can use same output path for Different Samples or Continue your job in Break Point.
        -t              threads number (default: 1)

    Additional parameters:
        -db             Add new sample to old database, old database path like /output/path/others/importDB_backup/
        
Written by Lu Liang
2025.7.11 at SHOU        
        ''')
        exit(0)   



def RunCMD(CMDin, info='', info_outfile=''):
    if info != '' and info_outfile != '':
        info_out = open(info_outfile, 'a+')
        info_out.write(info + '\t#Start\n')
        info_out.close()
    if CMDin.split()[0] == '#INFO':
        time_info = time.strftime("%Y-%m-%d %H:%M:%S INF: ", time.localtime())[2:]
        print(time_info + '\t'.join(CMDin.split()[1:]))
    else:
        time_info = time.strftime("%Y-%m-%d %H:%M:%S CMD: ", time.localtime())[2:]
        print(time_info + ' #Start ' + CMDin)
        os.system(CMDin)
        time_info = time.strftime("%Y-%m-%d %H:%M:%S CMD: ", time.localtime())[2:]
        print(time_info + ' # End  ' + CMDin)
    if info != '' and info_outfile != '':
        info_out = open(info_outfile, 'a+')
        info_out.write(info + '\t#End\n')
        info_out.close()
    

def BuildChrIndex(chr_name_list_in, ref_inp, out_idx_path):
    # Build index for each chromosome
    idx_log = ' > /dev/null 2>&1'
    def checkIDX(check_in_fa_path):
        # check if the index file exist
        if not os.path.exists(check_in_fa_path):
            return False
        request_idx_List = ['.amb', '.ann', '.bwt', '.fai', '.pac', '.sa']
        for ext_n in request_idx_List:
            if not os.path.exists(check_in_fa_path + i):
                return False
        if not os.path.exists(check_in_fa_path + '.dict'):
            return False
        return True
    
    g_seq_dict = {}
    seq_name = ''
    for line in open(ref_inp).readlines():
        if line[0] == '>':
            seq_name = line[1:].split()[0]
            g_seq_dict[seq_name] = [line]
        else:
            g_seq_dict[seq_name].append(line)
    
    idx_pool = Pool(int(t_num))
    if not checkIDX(out_idx_path + '/' + os.path.basename(ref_inp)):
        RunCMD('ln -s `realpath ' + ref_inp + '` ' + out_idx_path)
        raw_genome_path = out_idx_path + '/' + os.path.basename(ref_inp)
        idx_pool.apply_async(RunCMD, args=('samtools faidx ' + raw_genome_path + idx_log, ))
        idx_pool.apply_async(RunCMD, args=('bwa index ' + raw_genome_path + idx_log, ))
        idx_pool.apply_async(RunCMD, args=('gatk CreateSequenceDictionary -R ' + raw_genome_path + ' -O ' + os.path.splitext(raw_genome_path)[0] + '.dict' + idx_log, ))
    for n in chr_name_list_in:
        if not checkIDX(out_idx_path + '/' + n + '.fa'):
            outfile = open(out_idx_path + '/' + n + '.fa', 'w+')
            outfile.writelines(g_seq_dict[n])
            outfile.close()
            idx_pool.apply_async(RunCMD, args=('samtools faidx ' + out_idx_path + '/' + n + '.fa' + idx_log, ))
            idx_pool.apply_async(RunCMD, args=('bwa index ' + out_idx_path + '/' + n + '.fa' + idx_log, ))
            idx_pool.apply_async(RunCMD, args=('gatk CreateSequenceDictionary -R ' + out_idx_path + '/' + n + '.fa -O ' + out_idx_path + '/' + n + '.dict' + idx_log, ))
        # RunCMD('samtools faidx ' + referenceGenome + logout)
        # RunCMD('bwa index ' + referenceGenome + logout)
        # RunCMD('gatk CreateSequenceDictionary -R ' + referenceGenome + ' -O ' + os.path.splitext(referenceGenome)[0] + '.dict' + logout)
    idx_pool.close()
    idx_pool.join()
        
    
def SingleSampleGVCFProcess(runpath_in, ref_in, ref_seq_list_in, R1fq_in, R2fq_in, t_num_in, countID_in, run_db_path_in, backup_db_path_in):
    time_info_in = time.strftime("%Y-%m-%d %H:%M:%S Info: ", time.localtime())[2:]
    gatk_p_num = len(ref_seq_list_in)
    if gatk_p_num > 100:
        gatk_p_num = 100
    sample_start_time = time.time()
    RunCMD('mkdir -p ' + runpath_in)
    samplename = R1fq_in.split('/')[-1].split('_R1')[0]
    print(time_info_in, samplename, 'to gvcf Start', runpath_in)
    RunCMD('ln -s `realpath ' + R1fq_in + '` ' + runpath_in)
    RunCMD('ln -s `realpath ' + R2fq_in + '` ' + runpath_in)
    R1fq_in = runpath_in + '/' + R1fq_in.split('/')[-1]
    R2fq_in = runpath_in + '/' + R2fq_in.split('/')[-1]
    
    # 1. bwa alignment & sort
    logout = ' >> ' + runpath_in + '/run.log 2>&1'
    p_for_gunzip = Pool(2)
    rm_gunzip_list = []
    if '.gz' in R1fq_in:
        gun_name = R1fq_in.replace('.gz', '')
        gunzip_cmd = f'gunzip -c {R1fq_in} > {gun_name}'
        p_for_gunzip.apply_async(RunCMD, args=(gunzip_cmd, ))
        R1fq_in = gun_name
        rm_gunzip_list.append(R1fq_in)
    if '.gz' in R2fq_in:
        gun_name = R2fq_in.replace('.gz', '')
        gunzip_cmd = f'gunzip -c {R2fq_in} > {gun_name}'
        p_for_gunzip.apply_async(RunCMD, args=(gunzip_cmd, ))
        R2fq_in = gun_name
        rm_gunzip_list.append(R2fq_in)
    p_for_gunzip.close()
    p_for_gunzip.join()
    bwa_mem_shell = open(runpath_in + '/bwa_mem.sh', 'w+')
    bwa_mem_cmd = 'bwa mem -t ' + t_num_in + ' -R ' + '\'@RG\\tID:id' + str(countID_in) + '\\tPL:illumina\\tSM:' + samplename + '\' '+ ref_in + ' ' + R1fq_in + ' ' + R2fq_in + ' ' + '|samtools view -bS -F 4 -@ 10 -o ' + runpath_in + '/' + samplename + '.bam'
    bwa_mem_shell.write(bwa_mem_cmd)
    bwa_mem_shell.close()
    RunCMD('sh ' + runpath_in + '/bwa_mem.sh' + logout)
    # RunCMD('rm -rf ' + runpath_in + '/bwa_mem.sh')
    if len(rm_gunzip_list) > 0:
        for rm_gunzip in rm_gunzip_list:
            RunCMD('rm -rf ' + rm_gunzip)
    RunCMD('samtools sort -@ ' + t_num_in + ' -o ' + runpath_in + '/' + samplename + '.sorted.bam ' + runpath_in + '/' + samplename + '.bam' + logout)
    RunCMD('samtools index ' + runpath_in + '/' + samplename + '.sorted.bam' + logout)
    
    # Chromosomes in genome, split by seq then call
    # 2. samtools split bam by seq
    split_bam_outdir = runpath_in + '/splited_bam/'
    RunCMD('mkdir -p ' + split_bam_outdir)
    bam_split_pool = Pool(int(t_num_in))
    for seq in ref_seq_list_in:
        samtools_ext = 'samtools view -b ' + runpath_in + '/' + samplename + '.sorted.bam ' + seq + ' > ' + split_bam_outdir + seq + '.sorted.bam'
        bam_split_pool.apply_async(RunCMD, args=(samtools_ext, ))
    bam_split_pool.close()
    bam_split_pool.join()
    bam_idx_pool = Pool(int(t_num_in))
    for seq in ref_seq_list_in:
        samtools_idx = 'samtools index ' + split_bam_outdir + seq + '.sorted.bam'
        bam_idx_pool.apply_async(RunCMD, args=(samtools_idx, ))
    bam_idx_pool.close()
    bam_idx_pool.join()
    # RunCMD('rm -rf ' + runpath_in + '/' + samplename + '.sorted.bam ' + runpath_in + '/' + samplename + '.sorted.bam.bai ' + runpath_in + '/' + samplename + '.bam')
    
    # 4. mark PCR duplicates
    dup_t_num = str(int(int(t_num_in) / 10))
    P_for_MarkDup = Pool(10)
    for seq in ref_seq_list_in:
        markdupCMD = 'gatk MarkDuplicatesSpark ' + \
                        '-I ' + split_bam_outdir + seq + '.sorted.bam ' + \
                        '-O ' + split_bam_outdir + seq + '.sorted.markdup.bam ' + \
                        '-M ' + split_bam_outdir + seq + '.sorted.markdup_metrics.txt ' + \
                        '--spark-master local[' + dup_t_num + '] ' + logout
        P_for_MarkDup.apply_async(RunCMD, args=(markdupCMD, ))
    P_for_MarkDup.close()
    P_for_MarkDup.join()

    # for seq in ref_seq_list_in:
    #     RunCMD('rm -rf ' + split_bam_outdir + seq + '.sorted.bam')

    P_for_samtoolsIDX = Pool(int(t_num_in))
    for seq in ref_seq_list_in:
        dupidxCMD = 'samtools index ' + split_bam_outdir + seq + '.sorted.markdup.bam' + logout
        P_for_samtoolsIDX.apply_async(RunCMD, args=(dupidxCMD,))
        P_for_samtoolsIDX.apply_async(RunCMD, args=('rm -rf ' + split_bam_outdir + seq + '.sorted.bam ' + split_bam_outdir + seq + '.sorted.bam.bai', ))
    P_for_samtoolsIDX.close()
    P_for_samtoolsIDX.join()


    # 5. HaplotypeCaller processer
    gvcf_outdir = runpath_in + '/gvcf/'
    RunCMD('mkdir -p ' + gvcf_outdir)
    callSNP_pool = Pool(gatk_p_num)
    for seq in ref_seq_list_in:
        # ref_path =  os.path.dirname(ref_in) + '/' + seq + '.fa'
        callCMD = 'gatk --java-options \"-Xms12G -Xmx12G -XX:ParallelGCThreads=10\" HaplotypeCaller -ERC GVCF --native-pair-hmm-threads 10'
        callCMD += ' -I ' + split_bam_outdir + seq + '.sorted.markdup.bam'
        callCMD += ' -O ' + gvcf_outdir + seq + '.gvcf -R ' + ref_in + logout
        callSNP_pool.apply_async(RunCMD, args=(callCMD, ))
    callSNP_pool.close()
    callSNP_pool.join()
    # RunCMD('rm -rf ' + split_bam_outdir)

    # 6. GenomicsDBImport processer
    db_import_info = runpath_in + '/db_import_info/'
    RunCMD('rm -rf ' + db_import_info)
    RunCMD('mkdir -p ' + db_import_info)
    GenoDB_pool = Pool(gatk_p_num)
    while True:
        # line info: CMD    time    #Start/#End
        # check info file, if dbInport process running more than 3 hours, kill it and restart
        now_time = time.time()
        running_info = {}
        end_count = 0
        for file in os.listdir(db_import_info):
            for line in open(db_import_info + file).readlines():
                if '#Start' in line:
                    running_info[line.split('\t')[0]]['s'] = int(line.split('\t')[1])
                elif '#End' in line:
                    running_info[line.split('\t')[0]]['e'] = int(line.split('\t')[1])
                    end_count += 1
                elif '#Upload' in line:
                    running_info[line.split('\t')[0]] = {'u': int(line.split('\t')[1])}
        if end_count >= len(ref_seq_list_in):
            break
                
        for seq in ref_seq_list_in:
            db_path = run_db_path_in + '/' + seq + '_DB'
            GenoDBInfo = 'gatk GenomicsDBImport --genomicsdb-shared-posixfs-optimizations true --batch-size 50 ' + ' -L ' + seq + ' -V ' + gvcf_outdir + seq + '.gvcf --reader-threads 4' + logout
            info_outprint = GenoDBInfo + '\t' + str(int(now_time))
            info_outfile = db_import_info + seq + '_DB.log'
            if GenoDBInfo not in running_info:
                # run call cmd
                if os.path.exists(backup_db_path_in + '/' + seq + '_DB'):
                    GenoDBCMD = 'gatk GenomicsDBImport --genomicsdb-shared-posixfs-optimizations true --batch-size 50 --genomicsdb-update-workspace-path ' + db_path + ' -L ' + seq + ' -V ' + gvcf_outdir + seq + '.gvcf --reader-threads 4' + logout
                    info_out = open(info_outfile, 'w+')
                    info_out.write(GenoDBInfo + '\t' + str(int(now_time)) + '\t#Upload\n')
                    info_out.close()
                    GenoDB_pool.apply_async(RunCMD, args=(GenoDBCMD, info_outprint, info_outfile))
                else:
                    GenoDBCMD = 'gatk GenomicsDBImport --genomicsdb-shared-posixfs-optimizations true --batch-size 50 --genomicsdb-workspace-path ' + db_path + ' -L ' + seq + ' -V ' + gvcf_outdir + seq + '.gvcf --reader-threads 4' + logout
                    info_out = open(info_outfile, 'w+')
                    info_out.write(GenoDBInfo + '\t' + str(int(now_time)) + '\t#Upload\n')
                    info_out.close()
                    GenoDB_pool.apply_async(RunCMD, args=(GenoDBCMD, info_outprint, info_outfile))
            
            if GenoDBInfo in running_info and 's' in running_info[GenoDBInfo]:
                process_info = os.popen('ps -ef | grep ' + db_path + ' | grep gatk-package-4.3.0.0-local.jar|grep -v grep').readlines()
                if len(process_info) == 0:
                    info_out = open(info_outfile, 'a+')
                    info_out.write(GenoDBInfo + '\t' + str(int(now_time)) + '\t#End\n')
                    info_out.close()
                else:
                    if (now_time - running_info[GenoDBInfo]['s']) > 600:
                        RunCMD('kill ' + process_info[0].split()[1])
                        RunCMD('rm -rf ' + db_path)
                        if os.path.exists(backup_db_path_in + '/' + seq + '_DB'):
                            RunCMD('cp -r ' + backup_db_path_in + '/' + seq + '_DB ' + run_db_path_in)
                            GenoDBCMD = 'gatk GenomicsDBImport --genomicsdb-shared-posixfs-optimizations true --batch-size 50 --genomicsdb-update-workspace-path ' + db_path + ' -L ' + seq + ' -V ' + gvcf_outdir + seq + '.gvcf --reader-threads 4' + logout
                            info_out = open(info_outfile, 'w+')
                            info_out.write(GenoDBInfo + '\t' + str(int(now_time)) + '\t#Upload\n')
                            info_out.close()
                            GenoDB_pool.apply_async(RunCMD, args=(GenoDBCMD, info_outprint, info_outfile))
                        else:
                            GenoDBCMD = 'gatk GenomicsDBImport --genomicsdb-shared-posixfs-optimizations true --batch-size 50 --genomicsdb-workspace-path ' + db_path + ' -L ' + seq + ' -V ' + gvcf_outdir + seq + '.gvcf --reader-threads 4' + logout
                            info_out = open(info_outfile, 'w+')
                            info_out.write(GenoDBInfo + '\t' + str(int(now_time)) + '\t#Upload\n')
                            info_out.close()
                            GenoDB_pool.apply_async(RunCMD, args=(GenoDBCMD, info_outprint, info_outfile))
                                
        time.sleep(5)
    GenoDB_pool.close()
    GenoDB_pool.join()
    RunCMD('rm -rf ' + runpath_in)
    sample_end_time = time.time()
    # print(samplename, 'to gvcf Finish', runpath_in, 'Time:', str((sample_end_time - sample_start_time) / 3600)[:5], 'h')

    # 7. GenomeDB backup
    for seq in ref_seq_list_in:
        db_path = run_db_path_in + '/' + seq + '_DB'
        os.system('rm -rf ' + backup_db_path_in + '/' + seq + '_DB')
        os.system('cp -r ' + db_path + ' ' + backup_db_path_in)


def vcfFilter(raw_vcf_path_in, result_path_in, ref_seq_name_list_in, t_num_in, sample_count_in):
    logout = ' > ' + raw_vcf_path_in + '/filter_run.log 2>&1'
    # 5. get SNP and INDEL
    snpAndIndel_pool = Pool(int(t_num_in))
    for seq in ref_seq_name_list_in:
        if os.path.exists(raw_vcf_path_in + '/' + seq + '.vcf'):
            snpAndIndel_pool.apply_async(RunCMD, args=('gatk SelectVariants -V ' + raw_vcf_path_in + '/' + seq + '.vcf -O ' + raw_vcf_path_in + '/' + seq + '.snp.vcf --select-type-to-include SNP' + logout, ))
            snpAndIndel_pool.apply_async(RunCMD, args=('gatk SelectVariants -V ' + raw_vcf_path_in + '/' + seq + '.vcf -O ' + raw_vcf_path_in + '/' + seq + '.indel.vcf --select-type-to-include INDEL' + logout, ))
        pass
    snpAndIndel_pool.close()
    snpAndIndel_pool.join()

    # 6. filter vcf file
    filterVCF_pool = Pool(int(t_num_in))
    for seq in ref_seq_name_list_in:
        if os.path.exists(raw_vcf_path_in + '/' + seq + '.vcf'):
            filterVCF_pool.apply_async(RunCMD, args=('gatk VariantFiltration -O ' + raw_vcf_path_in + '/' + seq + '.snp.fil.vcf.temp -V ' + raw_vcf_path_in + '/' + seq + '.snp.vcf '
                + '--filter-expression \'QUAL < 30.0 || QD < 2.0 || FS > 60.0 ||  SOR > 4.0\' --filter-name lowQualFilter --cluster-window-size 10  --cluster-size 3 --missing-values-evaluate-as-failing' + logout, ))

            filterVCF_pool.apply_async(RunCMD, args=('gatk VariantFiltration -O ' + raw_vcf_path_in + '/' + seq + '.indel.fil.vcf.temp -V ' + raw_vcf_path_in + '/' + seq + '.indel.vcf '
                + '--filter-expression \'QUAL < 30.0 || QD < 2.0 || FS > 60.0 ||  SOR > 4.0\' --filter-name lowQualFilter --cluster-window-size 10  --cluster-size 3 --missing-values-evaluate-as-failing' + logout, ))
            pass
            
    filterVCF_pool.close()
    filterVCF_pool.join()

    # 7. filte PASS SNP/INDEL
    finalName1 = finalName + '_' + str(sample_count_in) + 'samples'
    RunCMD('grep \"#\" ' + raw_vcf_path_in + '/' + seq + '.snp.fil.vcf.temp >> ' + result_path_in + '/' + finalName1 + '.snp.fil.vcf')
    RunCMD('grep \"#\" ' + raw_vcf_path_in + '/' + seq + '.indel.fil.vcf.temp >> ' + result_path_in + '/' + finalName1 + '.indel.fil.vcf')
    for seq in ref_seq_name_list_in:
        if os.path.exists(raw_vcf_path_in + '/' + seq + '.vcf'):
            RunCMD('grep PASS ' + raw_vcf_path_in + '/' + seq + '.snp.fil.vcf.temp |grep -v \"#\" >> ' + result_path_in + '/' + finalName1 + '.snp.fil.vcf')
            RunCMD('grep PASS ' + raw_vcf_path_in + '/' + seq + '.indel.fil.vcf.temp |grep -v \"#\" >> ' + result_path_in + '/' + finalName1 + '.indel.fil.vcf')


def ReadOldDBSamplesName(old_db_json):
    json_dict = json.loads(open(old_db_json).read())
    old_db_sample_list = []
    for osp in json_dict['callsets']:
        old_db_sample_list.append(osp['sample_name'])
    return old_db_sample_list


if __name__ == '__main__':
    
    finalName = 'CallSNPRun_' + time.strftime("%Yy%mm%dd%Hh%Mm%Ss", time.localtime())[2:]
    if outputPath[-1] == '/':
        outputPath_name = outputPath[:-1].split('/')[-1]
    else:
        outputPath_name = outputPath.split('/')[-1]
    all_start_time = time.time()
    gatk_version = os.popen('gatk --version').read()
    if 'Genome Analysis Toolkit' not in gatk_version:
        print('### ERROR! gatk4 not install or not in PATH, please run \"conda activate gatk4\"')
        exit(0)
    if 'v4.' not in gatk_version:
        print('### ERROR! Please install GATK4.0 or later version')
        exit(0)
    # Read reference genome and get chromosome name
    ref_seq_name_list = []
    ref_seq_len_list = []
    all_genome_len = 0
    for line in open(referenceGenome).readlines():
        if line[0] == '>':
            ref_seq_len_list.append([line[1:].split()[0], 0])
        else:
            ref_seq_len_list[-1][1] += len(line.strip())
            all_genome_len += len(line.strip())
    need_len = 0
    for i in ref_seq_len_list:
        if i[1] / all_genome_len >= 0.001:
            ref_seq_name_list.append(i[0])
            need_len += i[1]
    if need_len / all_genome_len < 0.60:
        ref_seq_name_list = []
        for i in ref_seq_len_list:
            ref_seq_name_list.append(i[0])

    # Read fq files
    fqFileDict = {}
    for fqFile in os.listdir(fqFilePath):
        if fqFile.endswith('.fq') or fqFile.endswith('.fastq') or fqFile.endswith('.fq.gz') or fqFile.endswith('.fastq.gz'):
            if '_R1' in fqFile or '.R1' in fqFile:
                fq_name = fqFile.split('R1')[0]
                fq_path = fqFilePath + '/' + fqFile
            elif '_R2' in fqFile or '.R2' in fqFile:
                fq_name = fqFile.split('R2')[0]
                fq_path = fqFilePath + '/' + fqFile
            else:
                print('### ERROR fq name error', file)
                eixt(0)
            if fq_name[-1] == '_' or fq_name[-1] == '.':
                fq_name = fq_name[:-1]
            if fq_name not in fqFileDict:
                fqFileDict[fq_name] = [fq_path]
            else:
                fqFileDict[fq_name].append(fq_path)
    total_samples_count = len(fqFileDict)
    for i in fqFileDict:
        if len(fqFileDict[i]) != 2:
            print('### ERROR! Single fastq or repeat sample name', i, fqFileDict[i])
            exit(0)
    print(finalName, 'Start\n')
    
    resultPath = outputPath + '/result/'
    outputPath = outputPath + '/others/'
    # all_idx_path = outputPath + '/refDB/'
    backup_dbimport = outputPath + '/importDB_backup/'
    finished_sample_info = outputPath + '/finished_samples.txt'

    if not os.path.exists('/home/room/Lab_Users/gvcfDB_tmp/'):
        # create tmp gvcfDB
        os.system('mkdir -p /home/room/Lab_Users/gvcfDB_tmp/')
    else:
        # clean tmp gvcfDB
        for folder in os.listdir('/home/room/Lab_Users/gvcfDB_tmp/'):
            if len(os.popen('ps -ef|grep -v grep|grep gatk4_pipeline|grep ' + folder.split('_-')[1]).readlines()) == 0:
                os.system('rm -rf /home/room/Lab_Users/gvcfDB_tmp/' + folder)
        pass
    geno_DB_run = '/home/room/Lab_Users/gvcfDB_tmp/gvcfDB_-' + outputPath_name + '_-' + str(len(fqFileDict)) + '_samples_' + finalName.split('_')[1]
    print('Running CMD          : python3', ' '.join(sys.argv))
    print('Running DB path      :', geno_DB_run)
    # add new sample to old_db
    if old_db != '':
        print('Input database path   :', old_db)
        old_db_list = os.listdir(old_db)
        old_db_no_match_list = []
        for seq in ref_seq_name_list:
            if seq + '_DB' not in old_db_list:
                print('## ERROR ', old_db + '/' + seq + '_DB Not in reference genome ', referenceGenome )
                od_db_no_match_list.append(seq)
        for old_sp in ReadOldDBSamplesName(old_db + '/' + old_db_list[0] + '/callset.json'):
            if old_sp in fqFileDict:
                print('## ERROR Repeat sample name between Input database and input fastq:', old_sp)
                old_db_no_match_list.append(old_sp)
        if len(old_db_no_match_list) > 0:
            exit(0)
        else:
            RunCMD('cp -r ' + old_db + ' ' + outputPath)
    os.system('mkdir -p ' + geno_DB_run)
    finished_samples_list = []
    if os.path.exists(backup_dbimport):
        # check backup dbimport and del finished samples
        print('Recover Database......')
        if len(os.listdir(backup_dbimport)) == len(ref_seq_name_list):
            for seq in ref_seq_name_list:
                RunCMD('cp -r ' + backup_dbimport + '/' + seq + '_DB ' + geno_DB_run)
            new_sample_list = {}
            for line in open(finished_sample_info).readlines():
                finished_samples_list.append(line.split()[0])
            for sample in fqFileDict:
                if sample not in finished_samples_list:
                    new_sample_list[sample] = fqFileDict[sample]
            fqFileDict = new_sample_list
    else:
        os.system('mkdir -p ' + backup_dbimport)

    os.system('mkdir -p ' + resultPath)
    # RunCMD('mkdir -p ' + all_idx_path)
    fqFileDict_names = list(fqFileDict.keys())
    fqFileDict_names = sorted(fqFileDict_names)
    print('\nInput fq file include', len(fqFileDict), 'samples')
    for i in fqFileDict_names:
        print(i, fqFileDict[i])
    print('\nReference Genome in ', referenceGenome)
    print('Reference Genome include ', len(ref_seq_name_list), 'seqs')
    print('Threads          :', t_num)
    print('Output folder    :', outputPath)


    # 1.Build index
    # BuildChrIndex(ref_seq_name_list, referenceGenome, all_idx_path)
    # referenceGenome = all_idx_path + '/' + os.path.basename(referenceGenome)
    RunCMD('#INFO 0. Start build the index.')
    genomePath = os.path.dirname(referenceGenome)
    logout = ' > ' + genomePath + '/index_build.log'
    if genomePath == '':
        genomePath = '.'
    genomePathFileList = os.listdir(genomePath)
    requestList = ['.dict', '.amb', '.ann', '.bwt', '.fai', '.pac', '.sa']
    countRequest = 0
    if os.path.exists(os.path.splitext(referenceGenome)[0] + requestList[0]):
        countRequest += 1
    for n in requestList[1:]:
        if os.path.exists(referenceGenome + n):
            countRequest += 1
    if countRequest != len(requestList):
        RunCMD('samtools faidx ' + referenceGenome + logout)
        RunCMD('bwa index ' + referenceGenome + logout)
        RunCMD('gatk CreateSequenceDictionary -R ' + referenceGenome + ' -O ' + os.path.splitext(referenceGenome)[0] + '.dict' + logout)
        RunCMD('#INFO Build the index finished')
    else:
        RunCMD('#INFO Index already builded')


    # 2. Single sample gvcf
    if outputPath.split('/')[-1] == '':
        outputPath_name = outputPath.split('/')[-2]
    else:
        outputPath_name = outputPath.split('/')[-1]
    countID = len(finished_samples_list)
    running_ID = 0
    total_time = 0
    for sp in fqFileDict_names:
        start_time = time.time()
        countID += 1
        running_ID += 1
        op_in = outputPath + '/' + str(countID) + '_' + sp + '_' + str(len(fqFileDict) + len(finished_samples_list)) + 'total'
        if os.path.exists(op_in):
            RunCMD('rm -rf ' + op_in)
        R1_fq_name = ''
        R2_fq_name = ''
        for file in fqFileDict[sp]:
            if '_R1' in file or '.R1' in file:
                R1_fq_name = file
            elif '_R2' in file or '.R2' in file:
                R2_fq_name = file
            else:
                print('### ERROR fq name error', file)
                exit(0)
        # SingleSampleGVCFProcess(runpath_in, ref_in, ref_seq_list_in, R1fq_in, R2fq_in, t_num_in, countID_in, db_path_in):
        SingleSampleGVCFProcess(op_in, referenceGenome, ref_seq_name_list, R1_fq_name, R2_fq_name, t_num, countID, geno_DB_run, backup_dbimport)
        end_time = time.time()
        total_time = end_time - all_start_time
        time_info = time.strftime("%Y-%m-%d %H:%M:%S Info: ", time.localtime())[2:]
        remain_time = total_time * len(fqFileDict) / running_ID / 3600
        remain_time_d_h = str(remain_time / 24).split('.')[0] + ' d ' + str(float('0.' + str(remain_time / 24).split('.')[1]) * 24)[:3] + 'h'
        print(time_info, 
                countID, '/', total_samples_count, 
                'Single sample', sp, 'Finished', 
                'Time used:', str((end_time - start_time) / 3600)[:5], 'h,', 
                'Total Time:', str(total_time / 3600)[:5] + 'h', 
                'Estimated time remaining:', remain_time_d_h)
        finished_sample_info_out = open(finished_sample_info, 'a+')
        finished_sample_info_out.writelines(sp + '\t' + str((end_time - start_time) / 3600)[:5] + ' h\t' + time_info + '\n')
        finished_sample_info_out.close()
        
    # 3. joint call
    joint_out_path = outputPath + '/joint_call_vcf/'
    RunCMD('mkdir -p ' + joint_out_path)
    joint_call_pool = Pool(int(int(t_num) / 4))
    for seq in ref_seq_name_list:
        joint_call_pool.apply_async(RunCMD, args=('gatk --java-options \"-Xmx10G -XX:+UseParallelGC -XX:ParallelGCThreads=6\" GenotypeGVCFs -R ' + referenceGenome + ' -V gendb://' + outputPath + '/importDB_backup/' + seq + '_DB -O ' + joint_out_path + '/' + seq + '.vcf', ))
    joint_call_pool.close()
    joint_call_pool.join()
    end_time = time.time()
    total_time = end_time - all_start_time
    time_info = time.strftime("%Y-%m-%d %H:%M:%S Info: ", time.localtime())[2:]
    print('\n' + time_info + '###', len(ref_seq_name_list), 'sequences joint call finished, output to', joint_out_path, 'Total Time:', str(total_time / 3600)[:5], 'h\n')
    
    
    # 4. filter vcf
    vcfFilter(joint_out_path, resultPath, ref_seq_name_list, t_num, len(fqFileDict))
    end_time = time.time()
    total_time = end_time - all_start_time
    time_info = time.strftime("%Y-%m-%d %H:%M:%S Info: ", time.localtime())[2:]
    print('\n' + time_info + '### All finished, Total Time:', str(total_time / 3600)[:5], 'h\n')

import os
import argparse

parser = argparse.ArgumentParser(usage='python <PATH to this script>/runfsc.py [--fsc|F] [--dir|D] [--pre|P]',
                                 description='For running fsc 100 times')
parser.add_argument('--fsc', dest='f', help='directory of fsc installed', required=True)
parser.add_argument('--dir', dest='d', help='working directory', required=True)
parser.add_argument('--pre', dest='p', help='prefix of the file for which fsc need', required=True)
args = parser.parse_args()

os.chdir(f'{args.d}')

for i in range(1, 101):
    os.mkdir(f'./run{i}')
    os.system(f'cp ./{args.p}.tpl ./run{i}/')
    os.system(f'cp ./{args.p}.est ./run{i}/')
    os.system(f'cp ./{args.p}_jointMAFpop1_0.obs ./run{i}/')
    os.system(f'cp ./{args.p}_jointMAFpop2_0.obs ./run{i}/')
    os.system(f'cp ./{args.p}_jointMAFpop2_1.obs ./run{i}/')
    os.system(f'cp ./{args.p}_jointMAFpop3_0.obs ./run{i}/')
    os.system(f'cp ./{args.p}_jointMAFpop3_1.obs ./run{i}/')
    os.system(f'cp ./{args.p}_jointMAFpop3_2.obs ./run{i}/')
    os.chdir(f'./run{i}')
    os.system(f'{args.f}/fsc28 -t {args.p}.tpl -e {args.p}.est -m -0 -n 100000 -L 30 -s 0 -M -c 15 -B 15')
    os.chdir(f'{args.d}')
    if i == 100:
        #os.system('cp /mnt/disk2/Lab_Users/fangong/programs/selectbestrun.sh ./')
        os.system('nohup ./selectbestrun.sh > bestrun.log 2>&1 &')
print(f'''****************************************************
        Run {args.p} completed
****************************************************''')

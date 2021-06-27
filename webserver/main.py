# define imports
import argparse
import subprocess

# read inputs provided by user
parser = argparse.ArgumentParser()
parser.add_argument('--infile')
parser.add_argument('--NumReps', default=3)
parser.add_argument('--NumCond', default=4)
parser.add_argument('--cores', default=4)
parser.add_argument('--Experiment', default="ProtExample")
parser.add_argument('--paired', default="false")
parser.add_argument('--stat', default="true")
parser.add_argument('--secondcol', default="false")
parser.add_argument('--is_header', default="true")
parser.add_argument('--PreSetNumClustVSClust', default=0)
parser.add_argument('--PreSetNumClustStand', default=0)
parser.add_argument('--maxClust', default=20)
args = parser.parse_args()

f = open("properties.yml","a")
f.write(f'infile: {args.infile} \n')
f.write(f'NumReps: {args.NumReps} \n')
f.write(f'NumCond: {args.NumCond} \n')
f.write(f'cores: {args.cores} \n')
f.write(f'Experiment: {args.Experiment} \n')
f.write(f'paired: {args.paired} \n')
f.write(f'stat: {args.stat} \n')
f.write(f'secondcol: {args.secondcol} \n')
f.write(f'is_header: {args.is_header} \n')
f.write(f'PreSetNumClustVSClust: {args.PreSetNumClustVSClust} \n')
f.write(f'PreSetNumClustStand: {args.PreSetNumClustStand} \n')
f.write(f'maxClust: {args.maxClust} \n')
f.close()

fout = open("stdout.txt", "a")
ferr = open("stderr.txt", "a")
subprocess.run(["/srv/shiny-server/runVSClust.R", "properties.yml"], stdout=fout, stderr=ferr)

print("## VSClust - Results ")
print("View RPlots [here](Rplots.pdf) \n")
print(f"Download Stat output [here]({args.Experiment}statFileOut.csv) \n")
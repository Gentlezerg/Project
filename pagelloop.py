import pandas as pd
import subprocess

# import os   os.rename(filename,"indep"+filename)

profile = pd.read_csv("loopprofile.csv")

genelist = list(profile.columns.values)
genelist = genelist[1:len(genelist)]
length = len(genelist)
# length = 3
for i in range(length - 1):
    genei = genelist[i]
    for j in range(i + 1, length):
        genej = genelist[j]
        # print(genei,genej)
        filename = "pair"+genei+"and"+genej+".txt" # replace testpair
        logfile = filename+".log.txt" # replace testpair.txt.log.txt
        pairprofile = profile[["Unnamed: 0", genei, genej]]
        pairprofile.to_csv(filename, sep="	", header=False, index=False)
    
        # test # subprocess.call(["BayesTraitsV2","ClostRooted.trees","pairprofile.txt"],stdin=open("commandfile0.txt", 'r'))
        # subprocess.call(["BayesTraitsV2","ClostRooted.trees","testpair.txt"],stdin=open("commandfile0.txt", 'r'),shell=False)
        subprocess.call(["BayesTraitsV2_Quad", "ClostRooted.trees", filename], stdin=open("commandfile0.txt", 'r'),
                        shell=False)
        with open(logfile, 'r') as f:
            for line in f:
                pass
        indep = line.split()[1]
        # BayesTraitsV2_Quad
        subprocess.call(["BayesTraitsV2_Quad", "ClostRooted.trees", filename], stdin=open("commandfile1.txt", 'r'),
                        shell=False)
        with open(logfile, 'r') as f:
            for line in f:
                pass
        dep = line.split()[1]
    
        with open('pagelresult.txt', 'a') as results:
            results.write(",".join(map(str, [genei, genej, indep, dep])) + "\n")
        #print(genei, genej)



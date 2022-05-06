# DELETION QUANTIFIER V1.0

import csv
import os

def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)

basedir = "/home/USER/Desktop/WORKINGFOLDER"
setdirs = ["/home/USER/Desktop/WORKINGFOLDER/wt", "/home/USER/Desktop/WORKINGFOLDER/22", "/home/USER/Desktop/WORKINGFOLDER/co"]
basegene = "/home/USER/Desktop/WORKINGFOLDER/fragment.fasta"
datadir = basedir
resdir = os.path.join(basedir, "gapped_counts")

try:
    os.mkdir(resdir)
except:
    print("DIR EXISTS")

with open(basegene) as file:
    linesa = file.readlines()

subjecta = str(linesa[1])

filesa = os.listdir(setdirs[0])

files = []

for a in range(0, len(filesa)):
    if ".sam" in str(filesa[a]):
        files.append(str(filesa[a]))
for a in range(0, len(files)):

    with open(os.path.join(setdirs[0], files[a])) as file:
        lines = file.readlines()
    subject = [0]*(len(lines)-5)
    count = 0
    for b in range(5, len(lines)):
        subject[count] = lines[b]
        count = count + 1
    if a == 0:
        cigarfin = []
        lposfin = []
        flagfin = []
        seq1fin = []
        seq2fin = []
        readfin = []
        mapqfin = []
        basefin = []
    for b in range(0, len(subject)):
        for c in range(0, len(subject)):
            if c != b:
                inds = list(find_all(subject[b], '\t'))
                cigar1 = str(subject[b][inds[4]+1:inds[5]])
                lpos1 = str(subject[b][inds[2] + 1:inds[3]])
                flag1 = int(subject[b][inds[0] + 1:inds[1]])
                seq1 = str(subject[b][inds[8] + 1:inds[9]])
                read1 = str(subject[b][0:inds[0]])
                mapq1 = int(subject[b][inds[3]+1:inds[4]])

                inds = list(find_all(subject[c], '\t'))
                cigar2 = str(subject[c][inds[4]+1:inds[5]])
                lpos2 = str(subject[c][inds[2] + 1:inds[3]])
                flag2 = int(subject[c][inds[0] + 1:inds[1]])
                seq2 = str(subject[c][inds[8] + 1:inds[9]])
                read2 = str(subject[c][0:inds[0]])
                mapq2 = int(subject[c][inds[3]+1:inds[4]])

                readstring ="{}_|_{}".format(read1, read2)

                if readstring not in readfin:
                    if cigar1 == cigar2:
                        if "S" in cigar1 or "D" in cigar1: # s= skip, d=deletion
                            if lpos1 == lpos2:
                                lposfin.append(lpos1)
                                readfin.append(readstring)
                                cigarfin.append(cigar1)
                                flagfin.append("{}:{}".format(flag1,flag2))
                                seq1fin.append(seq1)
                                seq2fin.append(seq2)
                                readfin.append("{}_|_{}".format(read1, read2))
                                if mapq1 == mapq2:
                                    mapqfin.append(mapq1)
                                else:
                                    exit()
                                intlist = []
                                typelist = []
                                valuelist = []
                                totallength = 0

                                for d in range(0, len(cigar1)):
                                    try:
                                        int(cigar1[d])
                                        intlist.append(str(cigar1[d]))
                                    except:
                                        typelist.append(str(cigar1[d]))
                                        value = int(''.join(intlist))
                                        intlist = []
                                        valuelist.append(value)
                                        if str(cigar1[d]) != "S" or str(cigar1[d]) != "H":
                                            totallength = totallength + value
                                test = [0] * totallength
                                start = 0
                                for d in range(0, len(valuelist)):
                                    for e in range(0, valuelist[d]):
                                        subj = str(typelist[d])
                                        test[start] = subj
                                        start = start + 1
                                typelist = list.copy(test)
                                alnbase = subjecta[int(lpos1)-1:int(lpos1) + totallength]

                                base = list.copy(typelist)
                                sq1 = list.copy(typelist)
                                sq2 = list.copy(typelist)
                                countaln = 0
                                countsq1 = 0
                                countsq2 = 0
                                for d in range(0, len(typelist)):
                                    if typelist[d] == "S" or typelist[d] == "H":
                                        base[d] = "-"
                                        sq1[d] = seq1[countsq1]
                                        sq2[d] = seq2[countsq2]
                                        countsq1 = countsq1 + 1
                                        countsq2 = countsq2 + 1

                                    elif typelist[d] == "M":
                                        base[d] = alnbase[countaln]
                                        sq1[d] = seq1[countsq1]
                                        sq2[d] = seq2[countsq2]
                                        countaln = countaln + 1
                                        countsq1 = countsq1 + 1
                                        countsq2 = countsq2 + 1

                                    elif typelist[d] == "D":
                                        base[d] = alnbase[countaln]
                                        sq1[d] = " "
                                        sq2[d] = " "
                                        countaln = countaln + 1

                                    elif typelist[d] == "I":
                                        base[d] = " "
                                        sq1[d] = seq1[countsq1]
                                        sq2[d] = seq2[countsq2]
                                        countsq1 = countsq1 + 1
                                        countsq2 = countsq2 + 1

                                basefin.append(base)

rows = zip(readfin, cigarfin, lposfin, flagfin, seq1fin, seq2fin)

with open("/home/USER/Desktop/WORKINGFOLDER/wtstackDELS.csv", "w") as f:
    writer = csv.writer(f)
    for row in rows:
        writer.writerow(row)

with open(basegene) as file:
    linesa = file.readlines()
subjecta = str(linesa[1])

filesa = os.listdir(setdirs[1])

files = []
for a in range(0, len(filesa)):
    if ".sam" in str(filesa[a]):
        files.append(str(filesa[a]))
for a in range(0, len(files)):

    with open(os.path.join(setdirs[1], files[a])) as file:
        lines = file.readlines()
    subject = [0]*(len(lines)-5)
    count = 0
    for b in range(5, len(lines)):
        subject[count] = lines[b]
        count = count + 1
    if a == 0:
        cigarfin = []
        lposfin = []
        flagfin = []
        seq1fin = []
        seq2fin = []
        readfin = []
        mapqfin = []
        basefin = []
    for b in range(0, len(subject)):
        for c in range(0, len(subject)):
            if c != b:
                inds = list(find_all(subject[b], '\t'))
                cigar1 = str(subject[b][inds[4]+1:inds[5]])
                lpos1 = str(subject[b][inds[2] + 1:inds[3]])
                flag1 = int(subject[b][inds[0] + 1:inds[1]])
                seq1 = str(subject[b][inds[8] + 1:inds[9]])
                read1 = str(subject[b][0:inds[0]])
                mapq1 = int(subject[b][inds[3]+1:inds[4]])

                inds = list(find_all(subject[c], '\t'))
                cigar2 = str(subject[c][inds[4]+1:inds[5]])
                lpos2 = str(subject[c][inds[2] + 1:inds[3]])
                flag2 = int(subject[c][inds[0] + 1:inds[1]])
                seq2 = str(subject[c][inds[8] + 1:inds[9]])
                read2 = str(subject[c][0:inds[0]])
                mapq2 = int(subject[c][inds[3]+1:inds[4]])

                readstring ="{}_|_{}".format(read1, read2)

                if readstring not in readfin:
                    if cigar1 == cigar2:
                        if "S" in cigar1 or "D" in cigar1: # s= skip, d=deletion
                            if lpos1 == lpos2:
                                lposfin.append(lpos1)
                                readfin.append(readstring)
                                cigarfin.append(cigar1)
                                flagfin.append("{}:{}".format(flag1,flag2))
                                seq1fin.append(seq1)
                                seq2fin.append(seq2)
                                readfin.append("{}_|_{}".format(read1, read2))
                                if mapq1 == mapq2:
                                    mapqfin.append(mapq1)
                                else:
                                    exit()
                                intlist = []
                                typelist = []
                                valuelist = []
                                totallength = 0

                                for d in range(0, len(cigar1)):
                                    try:
                                        int(cigar1[d])
                                        intlist.append(str(cigar1[d]))
                                    except:
                                        typelist.append(str(cigar1[d]))
                                        value = int(''.join(intlist))
                                        intlist = []
                                        valuelist.append(value)
                                        if str(cigar1[d]) != "S" or str(cigar1[d]) != "H":
                                            totallength = totallength + value
                                test = [0] * totallength
                                start = 0
                                for d in range(0, len(valuelist)):
                                    for e in range(0, valuelist[d]):
                                        subj = str(typelist[d])
                                        test[start] = subj
                                        start = start + 1
                                typelist = list.copy(test)
                                alnbase = subjecta[int(lpos1)-1:int(lpos1) + totallength]

                                base = list.copy(typelist)
                                sq1 = list.copy(typelist)
                                sq2 = list.copy(typelist)
                                countaln = 0
                                countsq1 = 0
                                countsq2 = 0
                                for d in range(0, len(typelist)):
                                    if typelist[d] == "S" or typelist[d] == "H":
                                        base[d] = "-"
                                        sq1[d] = seq1[countsq1]
                                        sq2[d] = seq2[countsq2]
                                        countsq1 = countsq1 + 1
                                        countsq2 = countsq2 + 1

                                    elif typelist[d] == "M":
                                        base[d] = alnbase[countaln]
                                        sq1[d] = seq1[countsq1]
                                        sq2[d] = seq2[countsq2]
                                        countaln = countaln + 1
                                        countsq1 = countsq1 + 1
                                        countsq2 = countsq2 + 1

                                    elif typelist[d] == "D":
                                        base[d] = alnbase[countaln]
                                        sq1[d] = " "
                                        sq2[d] = " "
                                        countaln = countaln + 1

                                    elif typelist[d] == "I":
                                        base[d] = " "
                                        sq1[d] = seq1[countsq1]
                                        sq2[d] = seq2[countsq2]
                                        countsq1 = countsq1 + 1
                                        countsq2 = countsq2 + 1

                                basefin.append(base)

rows = zip(readfin, cigarfin, lposfin, flagfin, seq1fin, seq2fin)

with open("/home/USER/Desktop/WORKINGFOLDER/22stackDELS.csv", "w") as f:
    writer = csv.writer(f)
    for row in rows:
        writer.writerow(row)

with open(basegene) as file:
    linesa = file.readlines()
subjecta = str(linesa[1])

filesa = os.listdir(setdirs[2])

files = []
for a in range(0, len(filesa)):
    if ".sam" in str(filesa[a]):
        files.append(str(filesa[a]))
for a in range(0, len(files)):

    with open(os.path.join(setdirs[2], files[a])) as file:
        lines = file.readlines()
    subject = [0]*(len(lines)-5)
    count = 0
    for b in range(5, len(lines)):
        subject[count] = lines[b]
        count = count + 1
    if a == 0:
        cigarfin = []
        lposfin = []
        flagfin = []
        seq1fin = []
        seq2fin = []
        readfin = []
        mapqfin = []
        basefin = []
    for b in range(0, len(subject)):
        for c in range(0, len(subject)):
            if c != b:
                inds = list(find_all(subject[b], '\t'))
                cigar1 = str(subject[b][inds[4]+1:inds[5]])
                lpos1 = str(subject[b][inds[2] + 1:inds[3]])
                flag1 = int(subject[b][inds[0] + 1:inds[1]])
                seq1 = str(subject[b][inds[8] + 1:inds[9]])
                read1 = str(subject[b][0:inds[0]])
                mapq1 = int(subject[b][inds[3]+1:inds[4]])

                inds = list(find_all(subject[c], '\t'))
                cigar2 = str(subject[c][inds[4]+1:inds[5]])
                lpos2 = str(subject[c][inds[2] + 1:inds[3]])
                flag2 = int(subject[c][inds[0] + 1:inds[1]])
                seq2 = str(subject[c][inds[8] + 1:inds[9]])
                read2 = str(subject[c][0:inds[0]])
                mapq2 = int(subject[c][inds[3]+1:inds[4]])

                readstring ="{}_|_{}".format(read1, read2)

                if readstring not in readfin:
                    if cigar1 == cigar2:
                        if "S" in cigar1 or "D" in cigar1: # s= skip, d=deletion
                            if lpos1 == lpos2:
                                lposfin.append(lpos1)
                                readfin.append(readstring)
                                cigarfin.append(cigar1)
                                flagfin.append("{}:{}".format(flag1,flag2))
                                seq1fin.append(seq1)
                                seq2fin.append(seq2)
                                readfin.append("{}_|_{}".format(read1, read2))
                                if mapq1 == mapq2:
                                    mapqfin.append(mapq1)
                                else:
                                    exit()
                                intlist = []
                                typelist = []
                                valuelist = []
                                totallength = 0

                                for d in range(0, len(cigar1)):
                                    try:
                                        int(cigar1[d])
                                        intlist.append(str(cigar1[d]))
                                    except:
                                        typelist.append(str(cigar1[d]))
                                        value = int(''.join(intlist))
                                        intlist = []
                                        valuelist.append(value)
                                        if str(cigar1[d]) != "S" or str(cigar1[d]) != "H":
                                            totallength = totallength + value
                                test = [0] * totallength
                                start = 0
                                for d in range(0, len(valuelist)):
                                    for e in range(0, valuelist[d]):
                                        subj = str(typelist[d])
                                        test[start] = subj
                                        start = start + 1
                                typelist = list.copy(test)
                                alnbase = subjecta[int(lpos1)-1:int(lpos1) + totallength]

                                base = list.copy(typelist)
                                sq1 = list.copy(typelist)
                                sq2 = list.copy(typelist)
                                countaln = 0
                                countsq1 = 0
                                countsq2 = 0
                                for d in range(0, len(typelist)):
                                    if typelist[d] == "S" or typelist[d] == "H":
                                        base[d] = "-"
                                        sq1[d] = seq1[countsq1]
                                        sq2[d] = seq2[countsq2]
                                        countsq1 = countsq1 + 1
                                        countsq2 = countsq2 + 1

                                    elif typelist[d] == "M":
                                        base[d] = alnbase[countaln]
                                        sq1[d] = seq1[countsq1]
                                        sq2[d] = seq2[countsq2]
                                        countaln = countaln + 1
                                        countsq1 = countsq1 + 1
                                        countsq2 = countsq2 + 1

                                    elif typelist[d] == "D":
                                        base[d] = alnbase[countaln]
                                        sq1[d] = " "
                                        sq2[d] = " "
                                        countaln = countaln + 1

                                    elif typelist[d] == "I":
                                        base[d] = " "
                                        sq1[d] = seq1[countsq1]
                                        sq2[d] = seq2[countsq2]
                                        countsq1 = countsq1 + 1
                                        countsq2 = countsq2 + 1

                                basefin.append(base)

rows = zip(readfin, cigarfin, lposfin, flagfin, seq1fin, seq2fin)

with open("/home/USER/Desktop/WORKINGFOLDER/costackDELS.csv", "w") as f:
    writer = csv.writer(f)
    for row in rows:
        writer.writerow(row)

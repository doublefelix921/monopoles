import datetime

s=str(datetime.datetime.now())
print(s)
tstart=""

for n in [2,3,5,6,8,9]:
    tstart+=s[n]
tstart+="_"
for n in range(11,13):
    tstart+=s[n]
tstart+="."
for n in range(14,16):
    tstart+=s[n]
tstart+="."
for n in range(17,23):
    tstart+=s[n]

tstartsecs=0.0
tstartsecs+=float(s[0]+s[1]+s[2]+s[3])*31556926
tstartsecs+=float(s[5]+s[6])*2629743.83
tstartsecs+=float(s[8]+s[9])*86400
tstartsecs+=float(s[11]+s[12])*3600
tstartsecs+=float(s[14]+s[15])*60
tstartsecs+=float(s[17]+s[18])*1
tstartsecs+=float("0."+s[20]+s[21]+s[22]+s[23]+s[24])







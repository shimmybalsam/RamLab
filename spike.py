import re
file = open("check")
nonMnonconv = 0
nonMconv = 0
nonMwrong = 0
Mnonconv = 0
Mconv = 0
Mwrong = 0
all_lines = file.readlines()
primerR = re.compile("^(.)A(.)AT(.)AAATATTAATTAATTA")
# check1 = "CACATCAAATATTAATTAATTA"
# check2 = "aGACATCAAATATTAATTAATTA"
# check3 = "GACATCAAATATTAGTTAATTA"

for line in all_lines[1::4]:
    result1 = primerR.match(line[13:])
    if result1:
        crucial_places = [result1.group(1), result1.group(2), result1.group(3)]
        for group in crucial_places:
            if group == "C" or group == "T":
                if group == "C":
                    nonMnonconv+=1
                else:
                    nonMconv+=1
            else:
                nonMwrong+=1


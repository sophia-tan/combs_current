from out import out
opm = []
combs_db = []

for o in out:
    opm.append(o)

with open('combs_database.txt', 'r') as infile:
    for line in infile:
        combs_db.append(line.strip().split('.')[0][:4])

count = 0
for pdb in opm:
    if pdb in combs_db:
        count += 1
print(count)

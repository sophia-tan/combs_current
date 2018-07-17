opm = []
combs_db = []

with open('opm_entries.txt', 'r') as inF:
    for line in inF:
        opm.append(line.strip())

with open('combs_database.txt', 'r') as infile:
    for line in infile:
        combs_db.append(line.strip().split('.')[0][:4])

count = 0
for pdb in opm:
    if pdb in combs_db:
        count += 1
print(count)

        

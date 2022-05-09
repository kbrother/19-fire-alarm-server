with open("result.txt") as f:
    lines = f.read().split("\n")

for line in lines:
    if not line: continue
    check_time = line.split()[-1]
    print(check_time)

with open("memory.txt") as f:
    lines = f.read().split("\n")

for line in lines:
    if not line: continue
    memory = line.split()[-1]
    print(memory)

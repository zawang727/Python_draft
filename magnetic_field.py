import csv

def readMagnitCsv():
    with open('iris.csv', newline='') as csvfile:
        rows = csv.reader(csvfile, delimiter=':')
        print(rows)
        return rows
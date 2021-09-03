import csv

def readMagnitCsv(filepath):
    with open(filepath, encoding="utf-8-sig") as csvfile:
        rows = csv.reader(csvfile, delimiter=',')
        grid = []
        for row in rows:
            grid.append(row)
        #print(grid)
        return grid
    
readMagnitCsv('magnet_1_measurement_for_test.csv')
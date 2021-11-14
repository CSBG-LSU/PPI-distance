import sqlite3
con = sqlite3.connect('./example.db')
cur = con.cursor()

cur.execute('''DROP TABLE IF EXISTS stocks''')
cur.execute('''CREATE TABLE stocks
               (symbol text, price real)''')

cur.execute("INSERT INTO stocks VALUES ('RHAT', 35.14)")
cur.execute("INSERT INTO stocks VALUES ('RHAT', 35.14)")
cur.execute("INSERT INTO stocks VALUES ('RHAT', 35.14)")
cur.execute("INSERT INTO stocks VALUES ('ABB', 99.00)")
con.commit()


for row in cur.execute('SELECT * FROM stocks ORDER BY price'):
    print(row[0], row[1])
print('---------------------------------')
cur.execute('''DROP TABLE IF EXISTS temp_table''')
cur.execute("CREATE TABLE temp_table as SELECT DISTINCT * FROM stocks")
cur.execute("DELETE FROM stocks")
cur.execute("INSERT INTO stocks SELECT * FROM temp_table")
cur.execute('''DROP TABLE IF EXISTS temp_table''')
con.commit()

for row in cur.execute('SELECT * FROM stocks ORDER BY price'):
    print(row[0], row[1])

con.close()


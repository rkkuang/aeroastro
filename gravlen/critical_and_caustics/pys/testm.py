masses = []

Nrows = 4
for i in range(Nrows):
    # for j in range(Nrows):
        # masses += 
        # print([(-1)**(i+k) for k in range(Nrows)]) #mass

        minx = abs(Nrows//2)
        # print([k for k in range(-minx,Nrows-minx)]) #xs

        print([i-minx for k in range(Nrows)]) #ys

#THERE ARE MULTIPLE MATRICES INCLUDED IN THIS PROGRAM UNCOMMENT THE ONE YOU WISH TO USE#
import numpy as np
########################MATRIX TO SOLVE########################
A = np.array([[1.0, 1.0, 1.0, 1.0, 1.0],                      #
              [1.0, 1.0, 2.0, 3.0, 2.0],                      #
              [-1.0, 0.0, 2.0, 1.0, 1.0],                     #
              [3.0, 2.0, -1.0, 0.0, 1.0]])                    #
                                                              #
#A = np.array([[1.0, 2.0, 3.0, 4.0, 10.0],                     #
#              [-1.0, 1.0, 2.0, 3.0, 5.0],                     #
#              [1.0, -1.0, 1.0, 2.0, 3.0],                     #
#              [-1.0, 1.0, -1.0, 5.0, 4.0]])                   #
#                                                              #
#A = np.array([[1.0, -3.0, 7.0, 2.0],                          #
#              [2.0, 4.0, -3.0, -1.0],                         #
#              [-3.0, 7.0, 2.0, 3.0]])                         #
#                                                              #
#A = np.array([[3.0, -1.0, 3.0, 1.0, 6.0],                     #
#              [6.0, 0.0, 9.0, -2.0, 13.0],                    #
#              [-12.0, 0.0, -10.0, 5.0, -17.0],                #
#              [72.0, -8.0, 48.0, -19.0, 93.0]])               #
#                                                              #
#A = np.array([[3.0, 1.0, 4.0, -1.0, 7.0],                     #
#              [2.0, -2.0, -1.0, 2.0, 1.0],                    #
#              [5.0, 7.0, 14.0, -8.0, 20.0],                   #
#              [1.0, 3.0, 2.0, 4.0, -4.0]])                    #
###############################################################


######################################################################LOOP FOR GAUSSIAN##############################################################################
def Gauss(A):
    for p in range (n-1):                                       #Pivot Loop - changes the pivot position each time 
        if (userchoice == 1):                                   #If the user selected the 1st option
            print(end="")                                       #Nothing special happens, continue with normal gaussian
        elif (userchoice == 2):                                 #If the user selected the 2nd option
            partPiv(A, p)                                       #Calls Partial Pivot
        elif (userchoice == 3):                                 #If the user selected the 3rd option
            scaled(A, p)                                        #Calls Scaled Partial
        
        for r in range (p+1, n):                                #Changes row every time it loops
            if (A[(rv[p]),([p])] == 0):                         #If pivot is zero swap rows
                temp = rv[p]                                    #Saves the row number to move
                rv[p] = rv[r]                                   #Sets the next row as the row with a zero
                rv[p+1] = temp                                  #Puts the row with a zero in the next row
                continue                                        #Since the row below already has a zero go back to the begining of the row loop
           
            m = -A[(rv[r]),p]/A[(rv[p]),([p])]                  #Calculate m - (negative number to be zero) / (pivot)
            A[(rv[r]),p] = 0                                    #Sets the number below the pivot as zero
           
            for c in range (p+1, n+1):                          #Loop for columns
                A[(rv[r]),c] = A[(rv[r]),c] + m*A[rv[p],c]      #Multiplies the pivot row by m and adds it to the row below
    
    #######################################BACK SUBSTITUTION################################
    
    x = [0 for i in range(n)]                     #Creates x to store the solution

    for i in range(n-1, -1, -1):                  #Loop to access each row 
        x[i] = A[rv[i], n]/A[rv[i], i]            #Calculates the value of x as the last entry of the row over the numbers before it in the row
        
        for b in range(i-1, -1, -1):              #Loop to access each row and substitute in the calculated value of x
            A[rv[b], n] -= A[rv[b], i] * x[i]     #Modifies last column of the matrix based on the calculated values of x
    
    print()
    print("The solution to this matrix is: ", end='')
    for i in range(n):                     
        print("%.5f " % (x[i]), end='')           #Prints each number in the solution vector
    return
###################################################################################################################################################################

##############################################################DEFINES PARTIAL PIVOT################################################################################
def partPiv(A, p):
    maxEl = abs(A[rv[p],p])            #Takes the absolute value of the pivot
    maxRow = p                         #Sets the max row to the pivot
    for i in range (p+1, n):           #Loop to go through each number
        if (abs(A[rv[i],p]) > maxEl):  #Checks max of the column
            maxEl = abs(A[rv[i],p])    #Sets the max element to the current number
            maxRow = i                 #Sets the max row to the row the max element is in
        if (maxEl == 0):               #If the max element is zero return to thr top of the function
            continue                   
    if (maxRow > p):                   #Swaps row if nessecary
        tmp = rv[p]                    #Sets the current row to temp
        rv[p] = maxRow                 #Sets the current row to the max row
        rv[maxRow] = tmp               #Sets the row that the max row was originally in to the temp
    return                             #Returns the modified row vector
#####################################################################################################################################################################

################################################################SCALED PARTIAL PIVOT#################################################################################
def scaled(A, p):
    for i in range (p, n-1):                  #Searches through the columns
        maxR = 0                              #Sets the max of the row back to zero
        for j in range (n):                   #Searches through each row
            if (abs(A[rv[i], j]) > maxR):     #If the current entry is greater than the max it becomes the new max
                maxR = A[rv[i], j]
                si[rv[i]] = abs(maxR)         #Stores the max into a vector
                                              
    maxM = 0                                  #Sets the max m value back to zero
    maxI = 0                                  #Sets the index of m back to zero
    for i in range (p, n-1):                  #Goes through the rows
        m = (abs(A[rv[i], p])/si[rv[i]])      #Calculates m - (Entry in the matrix / The value S in si corresponding to the row)
        if(m > maxM):                         #If the calculated m is greater than the max then it becomes the new max
            maxM = m
            maxI = rv[i]                      #Stores the index of the max m
    if(maxI > p):
        tmp = rv[p]                  #Sets the current row to temp
        rv[p] = maxI                 #Sets the current row to the max row
        rv[maxI] = tmp               #Sets the row that the max row was originally in to the temp
    return rv                        #Returns the modified row vector
#####################################################################################################################################################################

#DECLARES VARIABLES AND VECTORS
n = len(A)                             #Gets length of the matrix
rv = [i for i in range (n)]            #Row vector used to swap rows
si = [0 for i in range (n)]            #Vector used to store the max values of a row
mi = [0 for i in range (n)]            #Vector used to store the max m in scaled partial pivot



print("Welcome to Project 2!")                                              #Prints Menu
print("")
print("This is the current Matrix being used in the code: \n",A)
print("")
print("Here is the menu options\n")
print(" ------Options Menu------\n"
          "1.Gaussian Elimination\n"
          "2.Partial Pivoting\n"
          "3.Scaled Partial Pivoting\n")


while True:                                                                 # Loops to check user-input    
    try:
        userchoice = int(input("Please select one by entering the assoiated numbers: ")) #Gets user input
    except ValueError:
        print("That is not a vaild, please try again.")                     # If the type is not a integer it will print an error text
                                                                            # better try again... Return to the start of the loop
        continue
    else:
        if (userchoice>3):                                                  #If userchoice is less greater than 3 then it will loop back to ask for input again
            continue
        elif (userchoice<1):                                                #If userchoice is less than 1 then it will loop back to ask for input again
            continue
        break                                                               # If type is right, the loop breaks
Gauss(A)
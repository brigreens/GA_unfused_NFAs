import csv


def listToString(s): 
    
    # initialize an empty string
    str1 = "," 
    
    # return string  
    return (str1.join(s))
        
        
list_csvs = ['acc_left_term_with_homo', 'acc_right_term_with_homo', 'acc_core_with_homo', 'don_core_with_homo', 'don_left_term_with_homo', 'don_right_term_with_homo']


for file in list_csvs:
    file_csv = file + '.csv'
    with open(file_csv, 'r') as csvfile:
        for line in csvfile:
            items = line.split(',')

            #new_line = []
            if len(items) == 3:
                if items[-1] == '\n':
                    list1 = items[:2]

            if len(items) == 6:
                if items[0] == '':
                    list2 = items[1:]
                    newline = list1 + list2
                    newline[-1] = newline[-1].strip()
                    print(newline)
                    newstr = listToString(newline)
                    print(newstr)
                    #append_list_as_row(r'C:\Users\bripe\OneDrive - University of Pittsburgh\GA\more_units\acc_core_new.csv', newline)
                    file_txt = file + '.txt'
                    with open(file_txt, "a+") as txtfile:
                        txtfile.write(newstr)
                        txtfile.write("\n")

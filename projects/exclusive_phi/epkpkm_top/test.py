my_list = (False, False, True, True)

def returnListResult(index_to_not_check, list_in):
    
    final_result=-1
    cc=0
    for ii in list_in:
        if cc != index_to_not_check:
            if( not ii ):
                return False
        cc+=1
    return True



outcome = returnListResult(1,my_list)
print(' out come is %s ' % (outcome))

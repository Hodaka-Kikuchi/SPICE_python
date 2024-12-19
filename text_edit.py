h_inc_entry = ces7.get()
    try:
        # 少数の場合
        if '/' not in h_inc_entry:
            h_inc = float(h_inc_entry)
        else:
            # 分数の場合
            h_inc = float(Fraction(h_inc_entry))
    except ValueError:
        pass
    
    k_inc_entry = ces8.get()
    try:
        # 少数の場合
        if '/' not in k_inc_entry:
            k_inc = float(k_inc_entry)
        else:
            # 分数の場合
            k_inc = float(Fraction(k_inc_entry))
    except ValueError:
        pass
    
    l_inc_entry = ces9.get()
    try:
        # 少数の場合
        if '/' not in l_inc_entry:
            l_inc = float(l_inc_entry)
        else:
            # 分数の場合
            l_inc = float(Fraction(l_inc_entry))
    except ValueError:
        pass
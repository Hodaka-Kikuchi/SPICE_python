cp_h_entry = cp_h.get()
    try:
        # 少数の場合
        if '/' not in cp_h_entry:
            cph = float(cp_h_entry)
        else:
            # 分数の場合
            cph = float(Fraction(cp_h_entry))
    except ValueError:
        pass
    
    cp_k_entry = cp_k.get()
    try:
        # 少数の場合
        if '/' not in cp_k_entry:
            cpk = float(cp_k_entry)
        else:
            # 分数の場合
            cpk = float(Fraction(cp_k_entry))
    except ValueError:
        pass
    
    cp_l_entry = cp_l.get()
    try:
        # 少数の場合
        if '/' not in cp_l_entry:
            cpl = float(cp_l_entry)
        else:
            # 分数の場合
            cpl = float(Fraction(cp_l_entry))
    except ValueError:
        pass
    
    cp = np.array([cph,cpk,cpl])
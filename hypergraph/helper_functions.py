def compute_L1(DE_dataframe, indices_from_pathway):
    temp = DE_dataframe.iloc[indices_from_pathway, 1:]

    up_reg = temp.where(temp > 0).sum(axis=0)
    down_reg = temp.where(temp < 0).sum(axis=0)
    total = temp.abs().sum(axis=0)

    sig_up_reg = temp.where(temp >= 1).sum(axis=0)
    sig_down_reg = temp.where(temp <= -1).sum(axis=0)
    sig_total = sig_up_reg + sig_down_reg.abs()

    # print(temp.where(temp >= 1))
    # print(temp.where(temp <= -1).count())  # count num genes sig down

    return (up_reg, down_reg, total, sig_up_reg, sig_down_reg, sig_total)


def dictionary_reverser(dictionary):
    new_dict = {}
    for key in dictionary.keys():
        for value in dictionary[key].split(" "):
            new_dict[value] = key
    return new_dict

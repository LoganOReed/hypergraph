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


def genes_to_reactions(network):
    # creates a dictionary where key is gene name
    # and values are lists of reactions that are effected
    # Note: some genes in network.genes have more than one name, we split
    # based on ' ' to have keys that correspond to exactly one name
    genes_to_reactions = {}
    for gene in network.genes:
        for name in gene.name.split(" "):
            genes_to_reactions[name] = gene.reaction.split(" ")
    return genes_to_reactions


def reactions_to_genes(self, network):
    reactions_to_genes = {}
    for gene in network.genes:
        for reaction in gene.reaction.split(" "):
            reactions_to_genes[reaction] = gene.name.split(" ")

    return reactions_to_genes


def get_well_ids_of_genes_in_network(self, clone_ORF_lookup, genes_to_reactions_dict):
        # in this section, we are looking up the well_ids to tie back to the
        # drug_data in order to get the effect on the reactions
        genes_in_network = []
        for key in genes_to_reactions_dict.keys():
            genes_in_network.append(key)

        well_id_in_network_list = []
        for gene in genes_in_network:
            try:
                # gene is prefixed by mtu:, clone_ORF_lookup excludes this prefix
                index = clone_ORF_lookup.iloc[:, 1].tolist().index(gene[4:])
                well_id_in_network_list.append(clone_ORF_lookup.iloc[index, 0])
            except ValueError:
                pass

        return well_id_in_network_list

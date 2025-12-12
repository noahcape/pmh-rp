from tlp import MultiLabeledTree


def test_small():
    edgelist = [
        (1, 3),
        (1, 13),
        (3, 11),
        (3, 14),
        (13, 12),
        (13, 6),
        (6, 4),
        (12, 5),
        (5, 10),
        (5, 7),
        (10, 2),
    ]

    labeling = {
        1: ["breast", "brain", "lung"],
        3: ["brain", "lung"],
        11: ["lung"],
        14: ["brain"],
        6: ["lung"],
        4: ["ribs"],
        5: ["liver"],
        10: ["liver"],
        7: ["kidney"],
        2: ["liver"],
    }

    clone_tree = MultiLabeledTree.from_edgelist(
        edgelist, labeling, ["breast", "brain", "lung", "liver", "ribs", "kidney"]
    )
    clone_tree.pull_down()
    character_set = list(set([item for value in labeling.values() for item in value]))

    for u,v in clone_tree.tree.edges:
        for c in character_set:
            for c2 in character_set:
                print(f"Migration from ({u}, {c}) to ({v}, {c2}) -> costs: {clone_tree.dist_f((u,v), c, c2)}")

    # print(clone_tree.dist_f((1, "1_breast"), "breast", "breast"))
    # print(clone_tree.dist_f((1, "1_breast"), "lung", "breast"))
    # print(clone_tree.dist_f((1, "1_breast"), "liver", "breast"))
    # print(clone_tree.dist_f((12, 5), "liver", "liver"))
    # print(clone_tree.dist_f((12, 5), "breast", "liver"))
    # print(clone_tree.leaf_f(2))
    # print(clone_tree.leaf_f(13))
    # clone_tree.solve_tlp()


def test_larger():
    edgelist = [
        ("grey", "red"),
        ("grey", "dark_brownC"),
        ("red", "mid_blue"),
        ("red", "light_green"),
        ("light_green", "dark_purple"),
        ("light_green", "dark_blue"),
        ("dark_purple", "orange"),
        ("dark_purple", "mid_green"),
        ("dark_purple", "dark_green"),
        ("orange", "pink"),
        ("orange", "light_purpleE"),
        ("mid_green", "goldF"),
        ("dark_green", "light_blue"),
        ("dark_green", "yellow"),
        ("dark_green", "light_brownI"),
    ]

    labeling = {
        "grey": ["seminal_vesicle", "prostate"],
        "red": ["left_pelvic_lymph_node_5", "seminal_vesicle	prostate"],
        "mid_blue": ["bladder", "right_pelvic_lymph_node_12", "seminal_vesicle"],
        "light_green": [
            "pelvic_lymph_node_7",
            "right_pelvic_lymph_node_12",
            "left_pelvic_lymph_node_8",
            "left_pelvic_lymph_node_5",
            "right_adrenal_gland",
            "bladder",
        ],
        "dark_purple": [
            "pelvic_lymph_node_7",
            "right_pelvic_lymph_node_12",
            "left_pelvic_lymph_node_8",
            "left_adrenal_gland",
            "left_pelvic_lymph_node_5",
            "right_adrenal_gland",
            "bladder",
        ],
        "dark_blue": ["bladder", "right_pelvic_lymph_node_12", "seminal_vesicle"],
        "orange": [
            "left_humerus_bone_marrow",
            "right_pelvic_lymph_node_12",
            "bladder",
            "left_adrenal_gland",
        ],
        "mid_green": ["right_adrenal_gland", "left_adrenal_gland"],
        "dark_green": ["left_pelvic_lymph_node_5", "left_pelvic_lymph_node_8"],
        "light_blue": ["left_pelvic_lymph_node_5", "pelvic_lymph_node_7"],
        "yellow": ["left_pelvic_lymph_node_5", "pelvic_lymph_node_7"],
        "pink": ["left_humerus_bone_marrow", "right_pelvic_lymph_node_12"],
        "dark_brownC": ["prostate"],
        "goldF": ["right_adrenal_gland"],
        "light_brownI": ["left_pelvic_lymph_node_8"],
        "light_purpleE": ["left_adrenal_gland"],
    }

    character_set = list(set([item for value in labeling.values() for item in value]))

    clone_tree = MultiLabeledTree.from_edgelist(edgelist, labeling, character_set)
    clone_tree.pull_down()

    print(clone_tree.dist_f(("grey", "1_prostate"), "seminal_vesicle", "prostate"))
    for label in character_set:
        print(f"{label} -> prostate\t", clone_tree.dist_f(("grey", "1_prostate"), f"{label}", "prostate"))


if __name__ == "__main__":
    test_small()

from __future__ import annotations


def constraints_indices_add(indices_list, num_add) -> list:
    """
    Add a specified number to each index in the given list of indices.

    Args:
        indices_list (list of int): A list of indices to be incremented.
        num_add (int): The number to add to each index in the list.

    Returns
        list of int: A new list with each index incremented by the specified number.
    """
    return [x + num_add for x in indices_list]

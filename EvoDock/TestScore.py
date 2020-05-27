

def GetScore(targetIndividual):
    """

    :return:
    """
    score = 0
    for aa in targetIndividual[0]:
        score += ord(aa)
    score = score / len(targetIndividual[0])

    return score


def GetScore2(targetIndividual):
    """

    :return:
    """
    score = 0
    if targetIndividual[0][0] == 'K':
        score += 10
    if targetIndividual[0][1] == 'I':
        score += 10
    if targetIndividual[0][2] == 'C':
        score += 10
    if targetIndividual[0][3] == 'N':
        score += 10
    if targetIndividual[0][4] == 'V':
        score += 10
    if targetIndividual[0][5] == 'N':
        score += 10

    return score

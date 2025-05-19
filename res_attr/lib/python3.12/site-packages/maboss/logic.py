"""utilitary functions to manipulate boolean expressions."""

import pyparsing as pp

def checkReserved(token):
    if token.lower() in ["true", "false", "not", "or", "and", "xor", "node"]:
        raise Exception("Name %s is reserved !" % token)
    return token

varName = (pp.Word(pp.alphas, pp.alphanums+'_'))
varName.setParseAction(lambda token: checkReserved(token[0]))

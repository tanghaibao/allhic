from hic/clm import newCLMFile, parse, parse_ids

var c = newCLMFile("test", "tests/test.clm")
echo c.parse_ids(true)
#c.parse()

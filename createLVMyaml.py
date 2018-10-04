verbos = True

dict = {"a": 1, "b": 2}
lvmdata = open("lvmdatsimTemplate.yml")
with open("lvmdatasimTemplace.yml", 'rb') as f:
    data = f.read()  # produces single string
for key in dict.keys():
    data.replace("__%s__placeholder"%key, dict[key])
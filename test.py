"""
Brief project testing framework
"""

class Test():
    def __init__(self, description, value, function, *args):
        self.function = function
        self.value = value
        self.description = description
        self.args = args
    def get_description(self):
        return f"Testing {self.function.__name__}({','.join([str(x) for x in self.args])}) == {self.value}"
    def check(self):
        return (self.function(*self.args) == self.value)

def my_add(a,b):
    return a+b

def test():
    t = Test("Test addition", 3, my_add, 1,2)
    print(t.get_description())
    print(t.check())
    s = Test("Expecte failure", 4, my_add, 1,1)
    print(s.check())

if __name__ == "__main__":
    test()

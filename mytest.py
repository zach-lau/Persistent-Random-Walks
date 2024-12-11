"""
Brief project testing framework
"""

class Test():
    def __init__(self, description, function, *args):
        self.function = function
        self.description = description
        self.args = args
    def get_description(self):
        return self.description
    def check(self):
        try:
            val = self.function(*self.args)
        except Exception as e:
            print("\033[31mFail\033[0m")
            print(e)
            return False
        if not val:
            print(f"\033[31mFailed with args {','.join([str(x) for x in self.args])}\033[0m")
        else:
            print("\033[32mPass\033[0m")
        return val

def my_add(a,b):
    return a+b

def test():
    t = Test("Test addition", lambda a,b: a == b, my_add(1,2), 3)
    print(t.get_description())
    print(t.check())
    s = Test("Expected failure", lambda a,b: a == b, my_add(1,1), 4)
    print(s.check())

if __name__ == "__main__":
    test()

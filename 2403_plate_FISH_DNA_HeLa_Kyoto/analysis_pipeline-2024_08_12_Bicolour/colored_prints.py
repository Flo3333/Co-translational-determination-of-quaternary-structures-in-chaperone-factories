class Color() :
    def __init__(self, text) -> None:
        self._ansi_code = '\033['
        self._ansi_code_end = '\033[00m'
        self.text = text

    def __str__(self) -> str:
        return self._ansi_code + self.text + self._ansi_code_end
    
    def __add__(self, other) -> str:
        if not hasattr(other, '__str__') : raise TypeError("Can only add to objects with str attributes")
        return self.__str__() + other.__str__()
    
    def __repr__(self) -> str:
        return self._ansi_code + self.text + self._ansi_code_end

class red(Color) :
    def __init__(self, text) -> None:
        super().__init__(text)

        self._ansi_code += '31m'

class green(Color) :
    def __init__(self, text) -> None:
        super().__init__(text)

        self._ansi_code += '32m'

class blue(Color) :
    def __init__(self, text) -> None:
        super().__init__(text)

        self._ansi_code += '34m'

class magenta(Color) :
    def __init__(self, text) -> None:
        super().__init__(text)

        self._ansi_code += '35m'

class gray(Color) :
    def __init__(self, text) -> None:
        super().__init__(text)

        self._ansi_code += '90m'

class cyan(Color) :
    def __init__(self, text) -> None:
        super().__init__(text)

        self._ansi_code += '96m'

class yellow(Color) :
    def __init__(self, text) -> None:
        super().__init__(text)

        self._ansi_code += '33m'
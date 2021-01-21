class DNAfrag:
    """
    This class is for DNA fragment connection.
    """
    def __init__(self,frag):
        """
        :param frag(int): initial input of DNA fragment
        """
        self.frag = frag
        self.head = frag[:15]
        self.tail = frag[-15:]
    def add_left(self,frag):
        """
        This function aims to
        :param frag:
        :return:
        """
        self.head = frag[:15]
        self.frag = frag+self.frag[15:]
    def add_right(self,frag):
        self.tail = frag[-15:]
        self.frag = self.frag+frag[15:]
    def check(self,frag):
        if self.head == frag[-15:]:
            self.add_left(frag)
            return True
        elif self.tail == frag[:15]:
            self.add_right(frag)
            return True
        else:
            return False
    def show_frag(self):
        return self.frag
    def show_head(self):
        return self.head
    def show_tail(self):
        return self.tail
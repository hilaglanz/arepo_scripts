class plotted_stream:
    pos_x = None
    pos_y = None
    vel_x = None
    vel_y = None
    def __init__(self, pos_x, pos_y, vel_x, vel_y):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.vel_x = vel_x
        self.vel_y = vel_y


class plotted_heatmap:
    pos_x = None
    pos_y = None
    slice = None
    def __init__(self, pos_x, pos_y, slice):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.slice = slice

class plotted_scatter:
    pos_x = None
    pos_y = None
    radius = None
    def __init__(self, pos_x, pos_y, radius):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.radius = radius



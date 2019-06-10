import cairo
from .interval import Interval

def draw(f):
    def wrapper(*args,**kwds):
        self = args[0]
        self.cr.move_to(-0.15/self.cr.cur_scale,0.00707/self.cr.cur_scale)
        info = "Line {}".format(self.cr.ct)
        self.cr.show_text(info)
        self.cr.ct += 1
        self.cr.move_to(0,0)
        if 'colour' in kwds.keys():
            self.cr.set_colour(kwds['colour'])
        out = f(*args,**kwds)
        self.cr.stroke()
        self.cr.vshift(0.04)
        self.cr.set_colour('black')
        return out
    return wrapper

def start_cairo(filename, height_ratio=1,scale=1):
    "Create a normalized pdf file CrContext with margins"
    sz = 500*scale
    left_margin = 80
    top_margin = 10
    right_margin = 10
    bottom_margin = 10

    hdim = sz+left_margin+right_margin
    vdim = sz*height_ratio+top_margin+bottom_margin

    ps = cairo.PDFSurface(filename, hdim, vdim)
    cr = CrContext(ps,scale=scale)
    cr.select_font_face("Inconsolata", cairo.FONT_SLANT_NORMAL,
                    cairo.FONT_WEIGHT_NORMAL)
    cr.set_font_size(0.015/scale)
    cr.translate(left_margin, top_margin)
    cr.set_line_width(1/sz)
    cr.scale(sz,sz)

    cr.set_colour('black')

    return cr

class Visual:
    def __init__(self,gn,*args,**kwds):
        self.cr = start_cairo(*args,**kwds)
        self.gn = gn

    @draw
    def net(self,iv_net,highlight=set(),**kwds):
        "Draw the net interval corresponding to alpha"
        for net_iv in iv_net.net:
            self.cr.interval(net_iv,label=self.gn.nb_set_type(net_iv),no_line=True,label_endpoints=True)
        self.cr.stroke()
        for iv in highlight:
            self.cr.iv_box(iv)
            self.cr.set_colour('red')
            self.cr.stroke()
            self.cr.set_colour('black')

    @draw
    def interval(self,iv_net,**kwds):
        "Draw the intervals of generation alpha"
        levels = [[]]
        for iv in iv_net.iv:
            done = False
            for level in levels:
                if level == [] or level[-1].b <= iv.a:
                    done = True
                    level.append(iv)
                    break
            if not done:
                levels.append([iv])

        for idx,level in enumerate(levels):
            for iv in level:
                self.cr.interval(iv)
            self.cr.vshift(0.02)

    @draw
    def nb_set(self):
        "Print the neighbour set types, and the associated index"
        for nb,idx in self.gn.nb_set_types.items():
            self.cr.show_text("{} : {}".format(idx,nb))
            self.cr.newline()
            self.cr.move_to(0,0)

    def show(self):
        "show the diagram on the page"
        self.cr.show_page()

class CrContext(cairo.Context):
    def __init__(self,*args,scale=1):
        self.ct = 0
        self.cur_scale = scale
        super(CrContext,self).__init__()
        self.font_dct = {'normal' : 0.015/self.cur_scale,
           'small' : 0.01/self.cur_scale,
           'tiny' : 0.005/self.cur_scale}
        self.colour_dct = {
                'black' : (0,0,0,1),
                'red' : (1,0,0,1),
                'gray' : (0.5,0.5,0.5,1),
                'blue' : (0,0,1,1),
                'green' : (0,1,0,1)
                }

    def set_font(self,val):
        try:
            self.set_font_size(self.font_dct[val])
        except KeyError:
            print(f"Invalid font size {val}")

    def set_colour(self,colour):
        try:
            self.set_source_rgba(*self.colour_dct[colour])
        except KeyError:
            print(f"Invalid colour {colour}")

    def newline(self):
        self.vshift(0.02)

    def interval(self, iv, label=None, no_line=False,label_endpoints=False):
        iv_height = 0.006/self.cur_scale
        translate = iv_height*self.cur_scale
        if label_endpoints:
            translate += 0.00707
        self.vshift(translate)
        mp = (iv.a+iv.b)/2
        if not no_line:
            self.move_to(iv.a,0)
            self.line_to(iv.b,0)
        self.move_to(iv.a,-iv_height)
        self.line_to(iv.a,iv_height)
        self.move_to(iv.b,-iv_height)
        self.line_to(iv.b,iv_height)

        if label is not None:
            self.move_to(mp,0)
            self.show_text_centred(label)

        if label_endpoints:
            for ep in (iv.a,iv.b):
                self.move_to(ep,-2*iv_height)
                self.set_font('tiny')
                self.show_text_centred(str(ep))
                self.set_font('normal')

        self.vshift(-translate)

    def iv_box(self, iv):
        iv_height = 0.006/self.cur_scale
        translate = iv_height*self.cur_scale
        self.vshift(translate)

        self.rectangle(iv.a,0,iv.b-iv.a,2*iv_height)

        self.vshift(-translate)

    def vshift(self,dist):
        "Normalized translation"
        self.translate(0,dist/self.cur_scale)

    def show_text_centred(self,text):
        ext = self.text_extents(text)
        x,y =  self.get_current_point()
        new_x = x-ext.width/2
        new_y = y+ext.height/2
        self.move_to(new_x,new_y)
        self.show_text(text)
        return (new_x,new_y,ext.width,ext.height)
        # x_bearing=0.0008799999999999988, y_bearing=-0.00624, width=0.002120000000000004, height=0.00624, x_advance=0.0049999999999999975, y_advance=0.0

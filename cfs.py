import math
from scipy.integrate import quad

n=0
axis=1
nfin = 10

def f(t):
    if 0 <= t < (math.pi*2/6):
        return [0, 6*t/(math.pi*2), 0, 0]
    elif (math.pi*2/6) <= t < (math.pi*2*2/6):
        return [0, 1, 6*(t-(math.pi*2/6))/(math.pi*2), 0]
    elif (math.pi*2*2/6) <= t < (math.pi*2*3/6):
        return [0, 1, 1, 6*(t-(math.pi*2*2/6))/(math.pi*2)]
    elif (math.pi*2*3/6) <= t < (math.pi*2*4/6):
        return [0, 1 - (6*(t-(math.pi*2*3/6)))/(math.pi*2), 1, 1]
    elif (math.pi*2*4/6) <= t < (math.pi*2*5/6):
        return [0, 0, 1 - (6*(t-(math.pi*2*4/6)))/(math.pi*2), 1]
    elif (math.pi*2*5/6) <= t < (math.pi*2*6/6):
        return [0, 0, 0, 1 - (6*(t-(math.pi*2*5/6)))/(math.pi*2)]
    else:
        if t<0:
            return f(t+math.pi*2)
        else:
            return f(t-math.pi*2)

def f_c(t, axis):
    return f(t)[axis]

def f_zero(t):
    global n, axis
    return f_c(t, axis)
def f_cos(t):
    global n, axis
    return f_c(t, axis)*math.cos(n*t)
def f_sin(t):
    global n, axis
    return f_c(t, axis)*math.sin(n*t)

zero_x = 0
cos_x = [0]
sin_x = [0]

zero_y = 0
cos_y = [0]
sin_y = [0]

zero_z = 0
cos_z = [0]
sin_z = [0]

for x in range(0, 10):
    axis = 1
    if x==0:
        res_a0, err_a0 = quad(f_zero, -math.pi, math.pi)
        a0 = res_a0 / (math.pi*2)
        zero_x = a0
        continue
    n = x
    res_cos, err_cos = quad(f_cos, -math.pi, math.pi)
    res_sin, err_sin = quad(f_sin, -math.pi, math.pi)
    an = res_cos / math.pi
    bn = res_sin / math.pi
    cos_x.append(an)
    sin_x.append(bn)

for y in range(0, 10):
    axis = 2
    if y==0:
        res_a0, err_a0 = quad(f_zero, -math.pi, math.pi)
        a0 = res_a0 / (math.pi*2)
        zero_y = a0
        continue
    n = y
    res_cos, err_cos = quad(f_cos, -math.pi, math.pi)
    res_sin, err_sin = quad(f_sin, -math.pi, math.pi)
    an = res_cos / math.pi
    bn = res_sin / math.pi
    cos_y.append(an)
    sin_y.append(bn)

for z in range(0, nfin):
    axis = 3
    if z==0:
        res_a0, err_a0 = quad(f_zero, -math.pi, math.pi)
        a0 = res_a0 / (math.pi*2)
        zero_z = a0
        continue
    n = z
    res_cos, err_cos = quad(f_cos, -math.pi, math.pi)
    res_sin, err_sin = quad(f_sin, -math.pi, math.pi)
    an = res_cos / math.pi
    bn = res_sin / math.pi
    cos_z.append(an)
    sin_z.append(bn)








class Quaternion:
    def __init__(self):
        pass

    def bar_quaternion(self, q):
        r, a, b, c = q
        return [r, -a, -b, -c]
    
    def quaternion_multiply(self, first, second):
        a, b, c, d = first
        e, f, g, h = second
        q = [
            (a*e - b*f - c*g - d*h),
            (a*f + b*e - c*h + d*g),
            (a*g + b*h + c*e - d*f),
            (a*h - b*g + c*f + d*e)
        ]
        return q
    def quaternion_vector_rotation(self, Cn, Dn, t):
        cos = math.cos(t/2)
        sin = math.sin(t/2)
        _, d, e, f = Dn
        Qn = [cos, sin*d, sin*e, sin*f]
        Qn_reverse = self.bar_quaternion(Qn)
        imp = self.quaternion_multiply(Qn, Cn)
        rotated = self.quaternion_multiply(imp, Qn_reverse)
        return rotated



class Ellipse(Quaternion):
    def __init__(self, a, b, c, d, e, f):
        Quaternion.__init__(self)
        
        self.data = [a, b, c, d, e, f]




        








        self.C = ( ( a**2 + c**2 + e**2 ) - ( b**2 + d**2 + f**2 ) ) / ( a*b + c*d + e*f ) # constant for calculating tp
        self.first_derivative_0_point = math.atan( (1/2) * ( -self.C + math.sqrt( self.C**2 + 4 ) ) )
        
        self.extreme_ts = []
        for n in range(5):
            value = self.first_derivative_0_point + n*(math.pi/2)
            if 0 <= value < 2*math.pi:
                self.extreme_ts.append( value )

        self.extreme_values = []
        for t in self.extreme_ts:
            self.extreme_values.append(self.ellipse(t))
        
        self.extreme_length = []
        for v in self.extreme_values:
            self.extreme_length.append(self.vector_length(v))
        
        self.major_t_index = 0 if self.extreme_length[0] > self.extreme_length[1] else 1
        self.minor_t_index = 1 - self.major_t_index

        self.major_t = self.extreme_ts[self.major_t_index]
        self.minor_t = self.extreme_ts[self.minor_t_index]

        self.major_radius = self.extreme_values[self.major_t_index]
        self.minor_radius = self.extreme_values[self.minor_t_index]
        
        self.major_length = self.extreme_length[self.major_t_index]
        self.minor_length = self.extreme_length[self.minor_t_index]

        self.pn_circle_length = (self.major_length + self.minor_length) / 2
        self.mn_circle_length = (self.major_length - self.minor_length) / 2

        self.pn_original = self.get_scaled_directional_vector(self.major_radius, self.pn_circle_length) # plus way rotation
        self.mn_original = self.get_scaled_directional_vector(self.major_radius, self.mn_circle_length) # minus way rotation
        self.t_original = self.major_t

        self.pDn = self.get_scaled_directional_vector(self.cross_product(self.extreme_values[1], self.extreme_values[0]), 1)
        self.mDn = self.get_scaled_directional_vector(self.pDn, -1)

        self.pCn = super().quaternion_vector_rotation(self.pn_original, self.pDn, -self.major_t)
        self.mCn = super().quaternion_vector_rotation(self.mn_original, self.mDn, -self.major_t)


        self.plus_way_rotation_cosine = self.pCn
        self.plus_way_rotation_sine = super().quaternion_vector_rotation(self.pCn, self.pDn, math.pi/2)

        self.minus_way_rotation_cosine = self.mCn
        self.minus_way_rotation_sine = super().quaternion_vector_rotation(self.mCn, self.mDn, math.pi/2)
        pass

    def ellipse(self, t):
        a, b, c, d, e, f = self.data
        X = a * math.cos(t) + b * math.sin(t)
        Y = c * math.cos(t) + d * math.sin(t)
        Z = e * math.cos(t) + f * math.sin(t)
        return [0, X, Y, Z]

    def vector_length(self, v):
        _, x, y, z = v
        return math.sqrt( x**2 + y**2 + z**2 )

    def get_scaled_directional_vector(self, v, sc):
        _, x, y, z = v
        scaling = sc / self.vector_length(v)
        return [0, scaling*x, scaling*y, scaling*z]
    
    def cross_product(self, u, v):
        _, u1, u2, u3 = u
        _, v1, v2, v3 = v
        cross = [0, u2*v3 - u3*v2, u3*v1 - u1*v3, u1*v2 - u2*v1]
        return cross

    def print_values(self):
        # print(f"Original : \nplus - {self.pn_original}\nminus - {self.mn_original}\nt - {self.t_original}")
        # print(f"plus-way rotation >>>>>>>>>>> \nCn : {self.pCn[1:]}\nDn : {self.pDn[1:]}\ncosine : {self.plus_way_rotation_cosine[1:]}\nsine : {self.plus_way_rotation_sine[1:]}")
        # print(f"minus-way rotation >>>>>>>>>>>>>> \nCn : {self.mCn[1:]}\nDn : {self.mDn[1:]}\ncosine : {self.minus_way_rotation_cosine[1:]}\nsine : {self.minus_way_rotation_sine[1:]}")
        print(f"plus-way rotation >>>>>>>>>>> \ncosine : {self.plus_way_rotation_cosine[1:]}\nsine : {self.plus_way_rotation_sine[1:]}")
        print(f"minus-way rotation >>>>>>>>>>>>>> \ncosine : {self.minus_way_rotation_cosine[1:]}\nsine : {self.minus_way_rotation_sine[1:]}")
        pass


    

# a, b, c, d, e, f = map(int, input().split())
# el = Ellipse(a, b, c, d, e, f)
# el.print_values()

el = []
for k in range(nfin):
    print("-------------------------------------------------------------------------------------------------------------"+str(k)+"---------------------------------")
    if k==0:
        el.append([zero_x, zero_y, zero_z])
        print(el[0])
        continue
    ell = Ellipse(cos_x[k], sin_x[k], cos_y[k], sin_y[k], cos_z[k], sin_z[k])
    el.append(ell)
    ell.print_values()

#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import pandas as pd


class MaterialProps:
    def __init__(self,props):
        # print(props.items())
        self.Ef    = props.e_f
        self.nu_f  = props.nu_f
        self.Em    = props.e_m
        self.nu_m  = props.nu_m
        self.vol_f = props.vol_f 

def get_material_properties():
    parser = argparse.ArgumentParser(description="Specify material properties")

    parser.add_argument("--e_f", type=float, default=241, help="Young's Modulus for the fiber (default: 241)")
    parser.add_argument("--nu_f", type=float, default=0.2, help="Poisson's Ratio for the fiber (default: 0.2)")
    parser.add_argument("--e_m", type=float, default=3.12, help="Young's Modulus for the matrix (default: 3.12)")
    parser.add_argument("--nu_m", type=float, default=0.38, help="Poisson's Ratio for the matrix (default: 0.38)")
    parser.add_argument("--vol_f", type=float, default=0.4, help="Volume Fraction of the fiber (default: 0.4)")
    parser.add_argument("--filename", type=str, default='pm_material_properties.csv', help="File name for saving props (default: pm_material_properties.csv)")

    args = parser.parse_args()

    return args

props = get_material_properties()
material_props = MaterialProps(props)
def print_values(material_props):
    print("User Input material properites")
    print(f"Ef:    {material_props.Ef}")
    print(f"nu_f:  {material_props.nu_f}")
    print(f"Em:    {material_props.Em}")
    print(f"nu_m:  {material_props.nu_m}")
    print(f"vol_f: {material_props.vol_f}")
    print("="*48)
print_values(material_props)
# defining parameters
e_f   = material_props.Ef
nu_f  = material_props.nu_f
e_m   = material_props.Em
nu_m  = material_props.nu_m
vol_f = material_props.vol_f


# In[2]:


#defining constants
s = {}
s[3] = 0.49247 - 0.47603 * vol_f - 0.02748 * vol_f * vol_f
s[6] = 0.36844 - 0.14944 * vol_f - 0.27152 * vol_f * vol_f
s[7] = 0.12346 - 0.32035 * vol_f + 0.23517 * vol_f * vol_f


# In[3]:


# calculating leme parameters
mu = {}
mu['f'] = e_f / (2 * (1 + nu_f))
mu['m'] = e_m / (2 * (1 + nu_m))

l_mbda = {}
l_mbda['m'] = e_m * nu_m / ((1 + nu_m) * (1 - 2 * nu_m))
l_mbda['f'] = e_f * nu_f / ((1 + nu_f) * (1 - 2 * nu_f))


# In[4]:


# calculating constants
a = mu['f'] - mu['m'] - 2 * mu['f'] * nu_m + 2 * mu['m'] * nu_f
b = - mu['m'] * nu_m + mu['f'] * nu_f + 2 * mu['m'] * nu_m * nu_f - 2 * mu['f'] * nu_m * nu_f
temp = (    mu['f'] - mu['m'] 
            + mu['f'] * nu_f 
            - mu['m'] * nu_m 
            + 2 * mu['m'] * nu_f 
            - 2 * mu['f'] * nu_m 
            + 2 * mu['m'] * nu_m * nu_f  
            - 2 * mu['f'] * nu_m * nu_f
            ) 
c = (mu['m'] - mu['f']) * temp
g = 2 - 2 * nu_m


# In[5]:


D = (
        a * s[3]**2 / (2 * mu['m']**2 * c)
    -   a * s[6] * s[3] / (mu['m']**2 * g * c)
    +   a * ( s[6]**2 - s[7]**2) / (2 * mu['m']**2 * g**2 * c) 
    +   s[3] * (b**2 - a**2) / (2 * mu['m'] * c**2)
    +   (s[6] * (a**2 - b**2) + s[7] * (a*b + b**2)) / (2 * mu['m'] * g * c**2)
    +   (a**3 - 2 * b**3 - 3 * a * b**2) / (8 * c**3)
)


# In[6]:


# calculating effective property Cstar
Cstar = {}

Cstar[44] = mu['m'] - vol_f / (- 2 * s[3]/mu['m'] + 1/(mu['m'] - mu['f']) + 4 * s[7]/(mu['m'] * (2 - 2 * nu_m)))

Cstar[66] = mu['m'] - vol_f / (- s[3]/mu['m'] + 1/(mu['m'] - mu['f']))


# In[7]:


temp = (
            s[3]**2 / mu['m']**2
            - 2 * s[6] * s[3] / ( mu['m']**2 * g)
            - a * s[3] / ( mu['m'] * c)
            + (s[6]**2 - s[7]**2) / ( mu['m']**2 * g**2 )
            + (a * s[6] + b * s[7]) / ( mu['m'] * g * c )
            + (a**2 - b**2) / ( 4 * c**2  )
        )
Cstar[11] = l_mbda['m'] + 2 * mu['m'] - vol_f /D * temp


# In[8]:


Cstar[12] = (
                l_mbda['m'] + 
                vol_f * b / D * ( 
                                    s[3] / ( 2 * c * mu['m']) 
                                    - (s[6] - s[7]) / ( 2 * c * mu['m'] * g) 
                                    - (a+b)/(4*c**2) 
                                    )
                )
Cstar[23] = (
                l_mbda['m'] + 
                vol_f  / D * ( 
                                    a * s[7] / ( 2 * g * c * mu['m']) 
                                    - ( b*a + b**2)/(4*c**2) 
                                    )
                )
Cstar[22] = (
                l_mbda['m'] + 
                2 * mu['m'] -
                vol_f  / D * ( 
                                    - a * s[3] / ( 2 * mu['m'] * c) 
                                    + a * s[6] / ( 2 * mu['m'] * g * c)  
                                    + (a**2 - b**2)/(4*c**2) 
                                    )
                )


# In[9]:


# calculating engineering constants
E1 = Cstar[11] - 2 * Cstar[12]**2/(Cstar[22] + Cstar[23])
E2 = ( 2 * Cstar[11] * Cstar[22] 
      + 2 * Cstar[11] * Cstar[23]
      - 4 * Cstar[12] * Cstar[12])
E2 *= (
        Cstar[22] - Cstar[23] + 2 * Cstar[44]
        )
E2 = E2 / (
            3 * Cstar[11] * Cstar[22]
            + Cstar[11] * Cstar[23]
            + 2 * Cstar[11] * Cstar[44]
            - 4 * Cstar[12] * Cstar[12]
            )
G12 = Cstar[66]
G13 = G12
nu12 = Cstar[12] / (Cstar[22] + Cstar[23])
nu13 = nu12
nu23 = (
            Cstar[11] * Cstar[22]
            + 3 * Cstar[11] * Cstar[23]
            - 2 * Cstar[11] * Cstar[44]
            - 4 * Cstar[12] * Cstar[12]
        ) 
nu23 = nu23 / (
                3 * Cstar[11] * Cstar[22] 
                +   Cstar[11] * Cstar[23]
                + 2 * Cstar[11] * Cstar[44]
                - 4 * Cstar[12] * Cstar[12]
                )
G23 = E2 / (2* (1 + nu23))


# In[10]:


results = {}
results['E1'] = E1
results['E2'] = E2
results['G12'] = G12
results['G13'] = G13
results['G23'] = G23
results['nu12'] = nu12
results['nu13'] = nu13
results['nu23'] = nu23


# In[11]:

print("="*48)
print("Output material properites")
for key, value in results.items():
    print(f'{key}: {value:.4f}')
print("="*48)


# In[ ]:


results = pd.DataFrame({
        'E1': [E1],
        'E2': [E2],
        'G12': [G12],
        'G13': [G13],
        'G23': [G23],
        'nu12': [nu12],
        'nu13': [nu13],
        'nu23': [nu23]
    })

    # Save the results to a CSV file
results.to_csv(props.filename, index=False)


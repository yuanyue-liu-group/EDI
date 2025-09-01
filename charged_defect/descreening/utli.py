import numpy as np
import io_xsf 
import lab

Bohr_R = 0.52917721067

def rho_expand(rho_array, target_fftn=np.array([241,241,241]), mode='constant'):
    pad_x = (target_fftn[0]-rho_array.shape[0])//2
    rest_x = (target_fftn[0]-rho_array.shape[0])%2

    pad_y = (target_fftn[1]-rho_array.shape[1])//2
    rest_y = (target_fftn[1]-rho_array.shape[1])%2

    print("pad_x: ", pad_x)
    print("rest_x: ", rest_x)
    print("pad_y: ", pad_y)
    print("rest_y: ", rest_y)

    extended_rho = np.pad(rho_array, 
                          pad_width=((pad_x, pad_x + rest_x), (pad_y, pad_y + rest_y), (0, 0)),
                            mode=mode, constant_values=0)
    print(rho_array.shape,'-->', extended_rho.shape)
    return extended_rho

def downsample_3d_array(arr, new_shape):
    s = [arr.shape[i] // new_shape[i] * new_shape[i] for i in range(3)]  # Find the nearest smaller shape divisible by new_shape
    print(s)
    arr = arr[:s[0], :s[1], :s[2]]  # Truncate the array to the new shape
    sh = (new_shape[0], arr.shape[0] // new_shape[0], 
          new_shape[1], arr.shape[1] // new_shape[1], 
          new_shape[2], arr.shape[2] // new_shape[2])  # Calculate the new shape
    print(sh)
    return arr.reshape(sh).mean(5).mean(1).mean(2)  # Take the mean along each dimension


def delta_v(arr1, arr2, align='None'):
    if align=='None':
        dv = arr2-arr1
    elif align=='vaccum':
        vacalign=np.average(arr2[:,:,0]-arr1[:,:,0])
        dv = arr2-arr1-vacalign
       
    return dv

def distance_array(nx, ny, nz, a,b,c, defect_loc='center'):
    x = np.linspace(0, a, nx)
    y = np.linspace(0, b, ny)
    z = np.linspace(0, c, nz)

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    if defect_loc == 'center':
        defect_loc = [0.5*a, 0.5*b, 0.5*c]

    return np.sqrt((X - defect_loc[0])**2 + (Y - defect_loc[1])**2 + (Z - defect_loc[2])**2)

def distance_array_2d(nx, ny, a,b, defect_loc='center'):
    x = np.linspace(0, a, nx)
    y = np.linspace(0, b, ny)

    X, Y = np.meshgrid(x, y, indexing='ij')
    
    if defect_loc == 'center':
        defect_loc = [0.5*a, 0.5*b]

    return np.sqrt(((X-0.5*Y) - (defect_loc[0]-0.5*defect_loc[1]))**2 + (np.sqrt(3)/2*Y - np.sqrt(3)/2*defect_loc[1])**2 )
    
def dat2npy(dat_inp, npy_outp):
    file = open(dat_inp)
    f = file.readlines()

    nx, ny, nz = int(f[1].split()[0]), int(f[1].split()[1]), int(f[1].split()[2])
    natom, ntype = int(f[1].split()[6]),int(f[1].split()[7])
    alat = float(f[2].split()[1])
    a1 = Bohr_R*alat*np.array([float(f[3].split()[0]),float(f[3].split()[1]),float(f[3].split()[2])])
    a2 = Bohr_R*alat*np.array([float(f[4].split()[0]),float(f[4].split()[1]),float(f[4].split()[2])])
    a3 = Bohr_R*alat*np.array([float(f[5].split()[0]),float(f[5].split()[1]),float(f[5].split()[2])])
    Nheader = 7+natom+ntype
    vl = ''.join(f[Nheader:])
    vn = [float(i) for i in vl.split()]
    __rho__ = np.array(vn).reshape([nz,ny,nx])
    rho = __rho__.transpose(2,1,0) 
    np.save(npy_outp, rho)
    





def read_dat(dat_inp):
    file = open(dat_inp)
    f = file.readlines()

    nx, ny, nz = int(f[1].split()[0]), int(f[1].split()[1]), int(f[1].split()[2])
    natom, ntype = int(f[1].split()[6]),int(f[1].split()[7])
    alat = float(f[2].split()[1])
    a1 = Bohr_R*alat*np.array([float(f[3].split()[0]),float(f[3].split()[1]),float(f[3].split()[2])])
    a2 = Bohr_R*alat*np.array([float(f[4].split()[0]),float(f[4].split()[1]),float(f[4].split()[2])])
    a3 = Bohr_R*alat*np.array([float(f[5].split()[0]),float(f[5].split()[1]),float(f[5].split()[2])])
    Nheader = 7+natom+ntype
    vl = ''.join(f[Nheader:])
    vn = [float(i) for i in vl.split()]
    __rho__ = np.array(vn).reshape([nz,ny,nx])
    rho = __rho__.transpose(2,1,0)
    
    return rho


def replace_rho_in_dat(ref_file, rho, out_file):
    """
    使用新的 rho 数组替换参考 .dat 文件中的电荷密度部分，保持 header 不变。
    """
    with open(ref_file, 'r') as f:
        lines = f.readlines()
    
    # 提取头部信息
    header_line = lines[1]
    nx, ny, nz = [int(i) for i in header_line.split()[:3]]
    expected_shape = (nx, ny, nz)

    if rho.shape != expected_shape:
        raise ValueError(f"rho shape {rho.shape} does not match expected shape {expected_shape} from {ref_file}")

    # 定位header结束位置
    natom = int(header_line.split()[6])
    ntype = int(header_line.split()[7])
    header_end = 7 + natom + ntype

    # 生成格式化电荷密度数据，每行写 4 个浮点数，空格分隔，换行结尾
    rho_flat = rho.transpose(2, 1, 0).flatten()
    data_lines = []
    for i, val in enumerate(rho_flat):
        data_lines.append(f"{val:17.10e}")
        if (i + 1) % 4 == 0:
            data_lines.append("\n")
        else:
            data_lines.append(" ")  # 添加空格以分隔数值

    # 确保最后一行也换行
    if (len(rho_flat) % 4) != 0:
        data_lines.append("\n")

    # 写入文件
    with open(out_file, 'w') as f_out:
        f_out.writelines(lines[:header_end])
        f_out.writelines(data_lines)

    print(f"Data has been written in: {out_file}")

if __name__ == "__main__": 
    cell = lab.mos2_tt
    old_sc = np.array([6,6,1])
    new_sc = np.array([12,12,1])

    Rho = io_xsf.Rho(cell['folder']+'rho_test.xsf')
    a1_ = list(map(float, Rho.a1))
    a2_ = list(map(float, Rho.a2))
    a1 = [i*2 for i in a1_]
    a2 = [i*2 for i in a2_]
    Rho.a1 = a1
    Rho.a2 = a2



    rho = Rho.rho

    rho_12x12 = rho_expand(rho, np.array([361,361,241]),mode='constant')


    rho_12x12_coarse = downsample_3d_array(rho_12x12, new_shape=np.array([180,180,120]))

    print(rho_12x12_coarse.shape)
    Rho.nx, Rho.ny, Rho.nz = rho_12x12_coarse.shape
    Rho.rho = rho_12x12_coarse

    io_xsf.write_xsf(Rho, xsf_type='rho', filedir=cell['folder']+'rho_6pad12_coarse.xsf')



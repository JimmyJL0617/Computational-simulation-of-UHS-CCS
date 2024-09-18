import os
import numpy as np
import pandas as pd
import cv2
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import interp1d as itp
from math import pi
from statistics import mean

Data_path = 'F://Server//'
first_dir = ['Zhenhuan']
Total_step = 70
res_threshold = 1000
max_gen, preser = [], []

def phase_sat(path,threshold):
    img = cv2.imread(path)
    gray = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
    H2 = img[:, :, 1]
    brine = img[:, :, 0]
    #Porosity
    pore_area = 0
    for k in range(len(img)):
        for q in range(len(img[0])):
            if not (H2[k][q] >= 60 and H2[k][q] <= 100) and not (brine[k][q] <= 136 and brine[k][q] >= 80):
                pore_area += 1
    porosity = pore_area/(len(H2)*len(H2[0]))
    #Total H2 and residual H2
    _, binary_H2 = cv2.threshold(H2, 110, 255, cv2.THRESH_BINARY)
    num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(binary_H2, 8, cv2.CV_32S)
    residual_area = 0
    H2_area = 0
    for j in range(1, num_labels):  # Start from 1 to ignore the background
        area = stats[j, cv2.CC_STAT_AREA]
        H2_area += area
        if area < threshold:
            residual_area += area
    res_H2_sat = residual_area/(len(H2)*len(H2[0])*porosity)
    return porosity, res_H2_sat

def run_process():
    for a in range(len(first_dir)):
        third_dir = ['LMH', 'HML', 'LHM', 'MLH', 'MHL', 'HLM']
        forth_dir = ['H2_aquifer_immis_compre_CA100_const', 'H2_aquifer_immis_compre_CA125_const', 'H2_aquifer_immis_compre_CA150_const']
        for c in range(len(third_dir)):
            for d in range(len(forth_dir)):
                path = Data_path + first_dir[a] + '//' + third_dir[c] + '//' + forth_dir[d] + '//'
                phase = []
                if os.path.exists(path):
                    imgs_H2 = img_prepro(path, third_dir[c])
                    count = len(imgs_H2)
                    for i in range(1, count+1):
                        if any('H2' in file and '.png' in file for file in os.listdir(path)):
                            porosity, res_H2_sat = phase_sat(path + imgs_H2[i-1], threshold=res_threshold)
                            if i<=10:
                                res_H2_sat = 0
                            df = pd.read_csv(path + 'Data_{}.csv'.format(round(i*0.001, 3)), names=['Time', 'X', 'Y', 'Z', 'p', 'alpha_H2', 'U1', 'U2', 'U3'], header=0)
                            H2_sat = df['alpha_H2'].mean()
                            Brine_sat = 1 - H2_sat
                            phase.append([i*1e-3, porosity, H2_sat, res_H2_sat, Brine_sat])
                            df_w = pd.DataFrame(phase, columns=['Time', 'Porosity', 'H2 saturation', 'Residual H2 saturation', 'Brine saturation'])
                            df_w.to_csv(Data_path + 'Zhenhuan//Plot//Processed_data' + '//' + third_dir[c] + '_' + forth_dir[d] + '.csv', index=False)
                            if i == count:
                                plot_data(df_w, title = first_dir[a] + '_' + third_dir[c] + '_' + forth_dir[d])
                                max_gen.append(df_w['Residual H2 saturation'].max())
                                preser.append(df_w['Residual H2 saturation'].iloc[-1])
                                df_w2 = pd.DataFrame({'Max generated residual H2': max_gen, 'Max preserved residual H2': preser})
                    pressure_ana(path, case_name = third_dir[c] + '_' + forth_dir[d], porosity=porosity, proc_file=Data_path + 'H2//Plot//Processed_data' + '//' + third_dir[c] + '_' + forth_dir[d] + '.csv', mode='space', col=c, row=d)
                        
    df_w2.to_csv(Data_path + 'Zhenhuan//Plot//' + 'Histagram' + '//' + 'H2 saturation.csv', index=False)
    get_max_conser('F://Server//Zhenhuan//Plot//Processed_data//')
    radial_histogram('F://Server//Zhenhuan//Plot//Histagram//')

def plot_data(file, title):
    #plt_file = pd.read_csv(Data_path + 'Processed_data//' + os.listdir(Data_path + 'Processed_data')[2])
    plt_file = file
    plt.figure(figsize=(12, 9))
    x = plt_file['Time']
    x = x*200
    y1 = plt_file['H2 saturation']
    y2 = plt_file['Residual H2 saturation']
    #y1_smooth = itp(x, y1, kind = 'nearest')(x)
    #y2_smooth = itp(x, y2, kind = 'nearest')(x)
    
    ax1 = plt.gca()
    ax1.tick_params(axis='both', labelsize=24)
    #ax1.scatter(x, y1)
    ax1.plot(x, y1, color='green', linewidth=5)
    
    for label in ax1.get_xticklabels() + ax1.get_yticklabels():
        label.set_fontweight('bold')
        label.set_fontname('Times New Roman')
        
    for label in ax1.get_yticklabels():
        label.set_color('green')
    
    ax2 = ax1.twinx()
    ax2.tick_params(axis='both', labelsize=24)
    #ax2.scatter(x, y2)
    ax2.plot(x, y2, color='blue', linewidth=5)
    
    for label in ax2.get_yticklabels():
        label.set_fontweight('bold')
        label.set_fontname('Times New Roman')
        label.set_color('blue')
    
    ax1.set_xlabel('Dimensionless time t*', fontdict={'name': 'Times New Roman', 'size': 28}, fontweight = 'bold', labelpad = 15)
    ax1.set_ylabel('Total H2', fontdict={'name': 'Times New Roman', 'size': 28}, color = 'green', fontweight = 'bold', labelpad = 15)
    ax1.set_ylim(-0.016, 0.8)
    ax2.set_ylabel('Residual trapped H2', fontdict={'name': 'Times New Roman', 'size': 28}, color = 'blue', fontweight = 'bold', labelpad = 15)
    ax2.set_ylim(-0.006, 0.3)
    ax1.legend()
    ax2.legend()
    plt.title(title)
    plt.show()
        
def get_xy(img, start_x=None, start_y=None, end_x=None, end_y=None):
    H2 = img[:, :, 1]
    #start_x, start_y, end_x, end_y = 0, 0, len(H2), len(H2[0])
    for a in range(len(H2)):
        if any (H2[a,:]!=87):
            start_y = a
            break
        else:
            continue
    for b in range(len(H2[0])):
        if any (H2[:,b]!=87):
            start_x = b
            break
        else:
            continue
    for c in range(start_y, len(H2)):
        if any (H2[c,:]!=87):
            continue
        else:
            end_y = c
            break
    for d in range(start_x, len(H2[0])):
        if any (H2[:,d]!=87):
            continue
        else:
            end_x = d
            break
    return start_x, start_y, end_x, end_y

def img_prepro(path, sec_dir):
    imgs_H2 = [file for file in os.listdir(path) if 'png' in file and 'H2' in file]
    #imgs_p = [file for file in os.listdir(path) if 'png' in file and not 'H2' in file]
        
    for i in range(len(imgs_H2)):
        im_H2 = cv2.imread(path + imgs_H2[i])
        #im_p = cv2.imread(path + imgs_p[i])
        
        if sec_dir=='LHM' or sec_dir=='HLM' or sec_dir=='HML' in path:
            im_flip_H2 = cv2.imread(path + imgs_H2[i])
            #im_flip_p = cv2.imread(path + imgs_p[i])
            im_flip_H2 = cv2.flip(im_flip_H2, 0)
            #im_flip_p = cv2.flip(im_flip_p, 0)
            cv2.imwrite(path + imgs_H2[i], im_flip_H2)
            #cv2.imwrite(path + imgs_p[i], im_flip_p)
        
        if len(im_H2) >= 370 and len(im_H2[0]) >= 130:
            start_x, start_y, end_x, end_y = get_xy(im_H2)
            cropped_im_H2 = im_H2[start_y:end_y, start_x:end_x]
            #cropped_im_p = im_p[start_y:end_y, start_x:end_x]
            cv2.imwrite(path + imgs_H2[i], cropped_im_H2)
            #cv2.imwrite(path + imgs_p[i], cropped_im_p)
                        
    return imgs_H2

#--------------------------------------------------------------------------------#

def pressure_ana(path, case_name,row,col,porosity,proc_file, mode='',area_num=45, time_split=2):
    
    csv_count = sum(1 for file in os.listdir(path) if file.endswith('.csv'))
    lists = [[] for _ in range(area_num)]
    p0_list = []

    for i in range(1, csv_count + 1):
        p = path + '//' + 'Data_{}.csv'.format(round(i*0.001, 3))
        df_p = pd.read_csv(p, names=['Time', 'X', 'Y', 'Z', 'p', 'alpha_H2', 'U1', 'U2', 'U3'], header=0)
        sorted_df = df_p.sort_values(by='X')
        if i % time_split == 0:
            p0 = sorted_df['p'].iloc[2]
            p0_list.append(p0)
                    
        seg_len = len(sorted_df) // area_num
        for k in range(area_num):
            start = k * seg_len
            end = (k+1) * seg_len if k < (area_num-1) else len(sorted_df)
            area_ave = sorted_df['p'][start:end].mean()
            '''
            if area_ave < 1.3e7:
                area_ave = area_ave - 1e7
            else:
                area_ave = area_ave - 1.5e7'''
            lists[k].append(area_ave)
            
    '''time = [(t+1)*0.02 for t in range(100)]
    data_out = {'Time': time}
    for z in range(area_num):
        data_out['A{}'.format(z)] = lists[z]
    df = pd.DataFrame(data_out) 
    df.to_csv('F://Server//H2//Hete//ver//Plot//Seg_pressure//' + case_name + '.csv', index=False)'''
    plt.figure(figsize = (12, 9))
    
    if mode=='time':
        time = [(t+1)*0.02 for t in range(70)]
        data_out = {'Time': time}
        for z in range(area_num):
            data_out['A{}'.format(z)] = lists[z]
        df = pd.DataFrame(data_out) 
        df.to_csv('F://Server//Zhenhuan//Hete//ver//Plot//Seg_pressure//' + case_name + '.csv', index=False)
        x = df['Time']
        y0 = (df['A0'] - df['A{}'.format(area_num/3-1)])/1000
        y1 = (df['A{}'.format(area_num/3)] - df['A{}'.format(2*area_num/3-1)])/1000
        y2 = (df['A{}'.format(2*area_num/3)] - df['A{}'.format(area_num-1)])/1000
        ax = plt.gca()
        ax.tick_params(labelsize=24)
        ax.plot(x, y0, x, y1, x, y2, linewidth=5)
        
        for label in ax.get_xticklabels():
            label.set_fontweight('bold')
            label.set_fontname('Times New Roman')
            
        for label in ax.get_yticklabels():
            label.set_fontweight('bold')
            label.set_fontname('Times New Roman')
        ax.set_xlabel('Dimensionless time', fontdict={'name': 'Times New Roman', 'size': 28}, fontweight = 'bold', labelpad = 15)
        ax.set_ylabel('ΔP (kPa)', fontdict={'name': 'Times New Roman', 'size': 28}, fontweight = 'bold', labelpad = 15)
        ax.set_ylim(-20,130)
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.legend(['1st layer', '2nd layer', '3rd layer'], prop={'family': 'Times New Roman', 'size': 24})
        plt.title(case_name)
        
    if mode=='space':
        time = [(t+1)*0.2 for t in range(70)]
        data_out = {'Time': time}
        for z in range(area_num):
            data_out['A{}'.format(z)] = lists[z]
        df = pd.DataFrame(data_out) 
        df.to_csv('F://Server//Zhenhuan//Seg_pressure//' + case_name + '.csv', index=False)
        x = []
        for t in range(area_num):
            x.append(t*(450/area_num) + (450/area_num)/2)
        
        sel_str = pd.read_csv('F://Server//Zhenhuan//Img_selection.csv', header=None).iloc[row+1, col+1]
        sel_spli = sel_str.split(',')
        y0 = df.iloc[4, 1:]-1e7
        y1 = df.iloc[eval(sel_spli[0])-1, 1:]-1e7
        y2 = df.iloc[eval(sel_spli[1])-1, 1:]-1e7
        y3 = df.iloc[eval(sel_spli[2])-1, 1:]-1e7
        y4 = df.iloc[eval(sel_spli[3])-1, 1:]-1e7  
        y5 = df.iloc[69, 1:]-1e7
            
        ax = plt.gca()
        ax.tick_params(labelsize=24)
        ax.plot(x, y0, x, y1, x, y2, x, y3, x, y4, x, y5, linewidth=5)
        
        for label in ax.get_xticklabels():
            label.set_fontweight('bold')
            label.set_fontname('Times New Roman')
            
        for label in ax.get_yticklabels():
            label.set_fontweight('bold')
            label.set_fontname('Times New Roman')
        ax.set_xlabel('Location (μm)', fontdict={'name': 'Times New Roman', 'size': 28}, fontweight = 'bold', labelpad = 15)
        ax.set_ylabel('ΔP (Pa)', fontdict={'name': 'Times New Roman', 'size': 28}, fontweight = 'bold', labelpad = 15)
        #ax.set_yscale('log')
        #ax.set_ylim(0,380)
        ax.xaxis.set_major_locator(MultipleLocator(50))
        #ax.yaxis.set_major_locator(MultipleLocator(40))
        ax.legend(['t* = 1', 't* = {}'.format(int(round(eval(sel_spli[0])/5, 0))), 't* = {}'.format(int(round(eval(sel_spli[1])/5, 0))), 't* = {}'.format(int(round(eval(sel_spli[2])/5, 0))), 't* = {}'.format(int(round(eval(sel_spli[3])/5, 0))), 't* = 14'], prop={'family': 'Times New Roman', 'size': 20})
        plt.title(case_name)
        
    if mode=='energy':
        df_H2_sat = pd.read_csv(proc_file)
        time = [(t*time_split)*0.02 for t in range(int(70/time_split))]
        data_out = {'Time': time}
        Ev = [z*100/(4.5*porosity*1e6) for z in p0_list]  #Mw/m3
        '''
        Ev = []
        for z in p0_list:
            if z < 1.3e7:
                z = z - 1e7
            else:
                z = z - 1.5e7
            Ev.append(z*100/(4.5*porosity*1e6))'''
        Em = []
        for s in range (len(p0_list)):
            if s<=5:
                Em.append(Ev[s]*1000/(9*df_H2_sat['H2 saturation'].iloc[s]))
            else:
                Em.append(Ev[s]*1000/(9*(df_H2_sat['H2 saturation'].iloc[s]-df_H2_sat['H2 saturation'].iloc[s-4:s].mean())))
        data_out['A1'] = Ev  #w/m^3
        data_out['A2'] = Em  #w/kgH2
        df = pd.DataFrame(data_out) 
        df.to_csv('F://Server//Zhenhuan//Energy_consum//' + case_name + '.csv', index=False)
        x = df['Time']
        y1 = df['A1']
        y2 = df['A2']
        ax1 = plt.gca()
        ax2 = ax1.twinx()
        ax1.tick_params(labelsize=24)
        ax2.tick_params(labelsize=24)
        ax1.plot(x, y1, color='orange', linewidth=5)
        ax2.plot(x, y2, color='green', linewidth=5)
        ax1.fill_between(x, y1, color='orange', alpha=0.1)
        ax2.fill_between(x, y2, color='green', alpha=0.1)
        
        for label in ax1.get_xticklabels():
            label.set_fontweight('bold')
            label.set_fontname('Times New Roman')
            
        for label in ax1.get_yticklabels():
            label.set_fontweight('bold')
            label.set_fontname('Times New Roman')
            label.set_color('orange')
            
        for label in ax2.get_yticklabels():
            label.set_fontweight('bold')
            label.set_fontname('Times New Roman')
            label.set_color('green')
            
        ax1.set_xlabel('Dimensionless time t*', fontdict={'name': 'Times New Roman', 'size': 28}, fontweight = 'bold', labelpad = 15)
        ax1.set_ylabel('Ev (Mw/m³)', fontdict={'name': 'Times New Roman', 'size': 28}, color = 'orange', fontweight = 'bold', labelpad = 15)
        ax2.set_ylabel('Em (kw/kgH₂)', fontdict={'name': 'Times New Roman', 'size': 28}, color = 'green', fontweight = 'bold', labelpad = 15)
        #ax2.set_yscale('log')
        #ax1.set_ylim(0,25)
        #ax2.set_ylim(8,500000)
        #ax.yaxis.set_major_locator(MultipleLocator(10))
        #ax.legend(['1st layer', '2nd layer', '3rd layer'], prop={'family': 'Times New Roman', 'size': 24})
        plt.title(case_name)
                           
    plt.show()

#--------------------------------------------------------------------------------#

def get_max_conser(path):
    probe = ['CA100', 'CA125', 'CA150']
    for t in range(len(probe)):
        max_conser = [[] for i in range(5)]
        Heter = []
        for file1 in os.listdir(path):
            if probe[t] in file1:
                df1 = pd.read_csv(path + file1)
                y_max = df1['Residual H2 saturation'].max()
                y_conser = df1['Residual H2 saturation'].iloc[-3:-1].mean()
                y_total = df1['H2 saturation'].iloc[-1]
                max_conser[0].append(y_max)
                max_conser[1].append(y_conser)
                max_conser[2].append(y_total)
                Heter.append(file1.split('_')[0] + '_' + file1.split('_')[-1].split('.')[0])
                
        for file2 in os.listdir('F://Server//Zhenhuan//Energy_consum//'):
            if probe[t] in file2:
                df2 = pd.read_csv('F://Server//Zhenhuan//Energy_consum//'+file2)
                y_em = df2['A1'].mean()
                y_ev = df2['A2'].mean()
                max_conser[3].append(y_em)
                max_conser[4].append(y_ev)
                
        data_out = {'Heterogeneity': Heter}
        data_out['Max'] = max_conser[0]
        data_out['conserved'] = max_conser[1]
        data_out['Total'] = max_conser[2]
        data_out['A1'] = max_conser[3]
        data_out['A2'] = max_conser[4]
        out = pd.DataFrame(data_out)
        out.to_csv('F://Server//Zhenhuan//Plot//Histagram//' + 'radial_histagram_{}.csv'.format(probe[t]), index=False)

def radial_histogram(csv_path):
    files = ['radial_histagram.csv']
    for csv_file in files:
        data = pd.read_csv(os.path.join(csv_path, csv_file))
        
        # Assuming 9 unique labels for 9 bins
        labels = data['Heterogeneity'].unique()
        num_labels = len(labels)
        angles = np.linspace(0, 2 * np.pi, num_labels, endpoint=False)
        
        # Setting up the plot
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 10))
        fig, ax1 = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 10))
        #fig, ax2 = plt.subplots(figsize=(16, 10))
        fig, ax2 = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 10))
        
        # Calculating width of each bar to Hver the entire bin
        width = 2 * pi / num_labels
        
        # Plot 'Max' data
        max_values = data['Max']
        ax.bar(angles, max_values, width=width, color='blue', alpha=0.6, edgecolor='black')
        
        # Plot 'conserved' data on top
        conserved_values = data['conserved']
        ax.bar(angles, conserved_values, width=width, color='red', alpha=0.8, edgecolor='black')
        
        # Adjust grid and labels
        ax.set_xticks(angles)  # Set ticks to the angles of the bars
        ax.set_xticklabels(labels, fontdict={'name': 'Times New Roman', 'size': 20}) 
        #ax.set_xticklabels([])  # Remove tick labels
        ax.set_yticklabels([])
        ax.tick_params(axis='y', which='major', labelsize=24)
        ax.yaxis.set_major_locator(MultipleLocator(0.05))
        ax.set_ylim(0,0.25)
        
        # Optional: Customize grid and axis for better visibility and alignment
        ax.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
        plt.title(csv_file)
        
        #ax1 = plt.gca()
        total_values = data['Total']
        ax1.bar(angles, total_values, width=width, color='red', alpha=0.6, edgecolor='black')
        # Adjust grid and labels
        ax1.set_xticks(angles)  # Set ticks to the angles of the bars
        ax1.set_xticklabels(labels, fontdict={'name': 'Times New Roman', 'size': 20})  # Remove tick labels
        ax1.set_yticklabels([])
        ax.tick_params(axis='y', which='major', labelsize=24)
        ax1.yaxis.set_major_locator(MultipleLocator(0.2))
        ax1.set_ylim(0,0.8)
        
        # Optional: Customize grid and axis for better visibility and alignment
        ax1.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
        
        plt.title(csv_file)
        '''
        #Plot Ev data
        max_values = data['A1']
        ax2.bar(angles, max_values, width=width, color='orange', alpha=0.6, edgecolor='black')
        
        # Plot Em data on top
        conserved_values = data['A2']
        ax2.bar(angles, conserved_values, width=width, color='green', alpha=0.8, edgecolor='black')
        
        # Adjust grid and labels
        ax2.set_xticks(angles)  # Set ticks to the angles of the bars
        ax2.set_xticklabels(labels, fontdict={'name': 'Times New Roman', 'size': 20}) 
        #ax.set_xticklabels([])  # Remove tick labels
        ax2.set_yticklabels([])
        ax2.tick_params(axis='y', which='major', labelsize=24)
        ax2.yaxis.set_major_locator(MultipleLocator(10000))
        #ax2.set_yscale('log')
        ax2.set_ylim(0,100000)
        
        # Optional: Customize grid and axis for better visibility and alignment
        ax2.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
        plt.title(csv_file)
        
        '''
        Ex = ['LMH_CA100_H2','HML_CA100_H2','LMH_CA125_H2','HML_CA125_H2','LMH_CA150_H2','HML_CA150_H2','LMH_CA100_CO2','HML_CA100_CO2','LMH_CA125_CO2','HML_CA125_CO2','LMH_CA150_CO2','HML_CA150_CO2']
        Ev_value = data['A1']
        Em_value = data['A2']
        
        ind = np.arange(len(Em_value))
        fig, ax2 = plt.subplots()
        
        # Bar plot for ax2
        ax2.bar(ind - 0.2, Ev_value, color='orange', width=0.4, alpha=0.85)
        ax2.set_ylabel('Ev (Mw/m³)', fontdict={'name': 'Times New Roman', 'size': 16}, color='orange', fontweight='bold', labelpad=15)
        ax2.tick_params(axis='y', labelsize=12, colors='orange')
        ax2.set_ylim(4, 16)
        
        # Create twin axis for the second bar plot
        ax3 = ax2.twinx()
        ax3.bar(ind + 0.2, Em_value, color='green', width=0.4, alpha=0.85)
        ax3.set_ylabel('Em (Mw/kgH₂)',fontdict={'name': 'Times New Roman', 'size': 16}, color='green', fontweight='bold', labelpad=15)
        ax3.tick_params(axis='y', labelsize=12, colors='green')
        ax3.set_yscale('log')
        ax3.set_ylim(100, 25000)
        
        # Set x-axis and labels
        ax2.set_xticks(ind)
        ax2.set_xticklabels(Ex, fontname='Times New Roman', fontsize=10, fontweight='bold', rotation=45)
        #plt.title('Your Title Here')  # Replace with your actual title
    plt.show()

#--------------------------------------------------------------------------------#

run_process()   
#get_max_conser('F://Server//H2//Plot//Processed_data//')
#radial_histogram('F://Server//H2//Plot//Histagram//')
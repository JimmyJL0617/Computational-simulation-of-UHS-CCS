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

Data_path = 'F://Server//CO2//'
first_dir = ['Hete']#, 'Control']
Total_step = 100
res_threshold = 1000
max_gen, preser = [], []

def phase_sat(path,threshold):
    img = cv2.imread(path)
    gray = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
    CO2 = img[:, :, 2]
    brine = img[:, :, 0]
    #Porosity
    pore_area = 0
    for k in range(len(img)):
        for q in range(len(img[0])):
            if not (CO2[k][q] >= 50 and CO2[k][q] <= 100) and not (brine[k][q] <= 136 and brine[k][q] >= 80):
                pore_area += 1
    porosity = pore_area/(len(CO2)*len(CO2[0]))
    #Total CO2 and residual CO2
    _, binary_CO2 = cv2.threshold(CO2, 83, 255, cv2.THRESH_BINARY)
    num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(binary_CO2, 8, cv2.CV_32S)
    residual_area = 0
    CO2_area = 0
    for j in range(1, num_labels):  # Start from 1 to ignore the background
        area = stats[j, cv2.CC_STAT_AREA]
        CO2_area += area
        if area < threshold:
            residual_area += area
    res_CO2_sat = residual_area/(len(CO2)*len(CO2[0])*porosity)
    return porosity, res_CO2_sat

def run_process():
    for a in range(len(first_dir)):
        if a == 0:
            second_dir = ['ver']#, 'hor']
            third_dir = ['LMH','HML']# 'LHM', 'MHL', 'MLH', 'HML', 'HLM']
            forth_dir = ['CO2_aquifer_immisible_CA100g', 'CO2_aquifer_immisible_CA125g', 'CO2_aquifer_immisible_CA150g']
            for b in range(len(second_dir)):
                for c in range(len(third_dir)):
                    for d in range(len(forth_dir)):
                        path = Data_path + first_dir[a] + '//' + second_dir[b] + '//' + third_dir[c] + '//' + forth_dir[d] + '//'
                        phase = []
                        if os.path.exists(path):
                            imgs_CO2, imgs_p = img_prepro(path, third_dir[c])
                            count = len(imgs_CO2)
                            for i in range(1, count+1):
                                if any('CO2' in file and '.png' in file for file in os.listdir(path)):
                                    porosity, res_CO2_sat = phase_sat(path + imgs_CO2[i-1], threshold=res_threshold)
                                    if i<=10:
                                        res_CO2_sat = 0
                                    df = pd.read_csv(path + 'Data_{}.csv'.format(round(i*0.0001, 4) if i != 0 else 0))
                                    CO2_sat = df['alpha.CO2'].mean()
                                    Brine_sat = 1 - CO2_sat
                                    phase.append([i*1e-4, porosity, CO2_sat, res_CO2_sat, Brine_sat])
                                    df_w = pd.DataFrame(phase, columns=['Time', 'Porosity', 'CO2 saturation', 'Residual CO2 saturation', 'Brine saturation'])
                                    df_w.to_csv(Data_path + 'Hete//ver//Plot//Processed_data' + '//' + third_dir[c] + '_' + forth_dir[d] + '.csv', index=False)
                                    if i == count:
                                        plot_data(df_w, title = first_dir[a] + '_' + second_dir[b] + '_' + third_dir[c] + '_' + forth_dir[d])
                                        max_gen.append(df_w['Residual CO2 saturation'].max())
                                        preser.append(df_w['Residual CO2 saturation'].iloc[-1])
                                        df_w2 = pd.DataFrame({'Max generated residual CO2': max_gen, 'Max preserved residual CO2': preser})
                            pressure_ana(path, case_name = third_dir[c] + '_' + forth_dir[d], porosity=porosity, proc_file=Data_path + 'Hete//ver//Plot//Processed_data' + '//' + third_dir[c] + '_' + forth_dir[d] + '.csv', mode='space', col=c, row=d)
                        else:
                            continue
        if a == 1:
            second_dir = ['Low', 'Med', 'High']
            third_dir = ['CO2_aquifer_immisible_CA100', 'CO2_aquifer_immisible_CA125', 'CO2_aquifer_immisible_CA150']
            for b in range(len(second_dir)):
                for c in range(len(third_dir)):
                    path = Data_path + first_dir[a] + '//' + second_dir[b] + '//' + third_dir[c] + '//'
                    phase = []
                    if os.path.exists(path):
                        imgs_CO2, imgs_p = img_prepro(path, second_dir[b])
                        count = len(imgs_CO2)
                        for i in range(1, count+1):
                            if any('CO2' in file and '.png' in file for file in os.listdir(path)):
                                porosity, res_CO2_sat = phase_sat(path + imgs_CO2[i-1], threshold=res_threshold)
                                if i<=10:
                                    res_CO2_sat = 0
                                df = pd.read_csv(path + 'Data_{}.csv'.format(round(i*0.0001, 4) if i != 0 else 0))
                                CO2_sat = df['alpha.CO2'].mean()
                                Brine_sat = 1 - CO2_sat
                                phase.append([i*1e-4, porosity, CO2_sat, res_CO2_sat, Brine_sat])
                                df_w = pd.DataFrame(phase, columns=['Time', 'Porosity', 'CO2 saturation', 'Residual CO2 saturation', 'Brine saturation'])
                                df_w.to_csv(Data_path + 'Hete//ver//Plot//Processed_data' + '//' + third_dir[c] + '.csv', index=False)
                                if i == count:
                                    plot_data(df_w, title = first_dir[a] + '_' + second_dir[b] + '_' + third_dir[c])
                                    max_gen.append(df_w['Residual CO2 saturation'].max())
                                    preser.append(df_w['Residual CO2 saturation'].iloc[-1])
                                    df_w2 = pd.DataFrame({'Max generated residual CO2': max_gen, 'Max preserved residual CO2': preser})
                        pressure_ana(path, case_name = second_dir[b] + '_' + third_dir[c], porosity=porosity, proc_file=Data_path + 'Hete//ver//Plot//Processed_data' + '//' + third_dir[c] + '.csv', mode='space', col=b, row=c)
                    else:
                        continue
    df_w2.to_csv(Data_path + 'Hete//ver//Plot//' + 'Histagram' + '//' + 'CO2 saturation.csv', index=False)
    get_max_conser('F://Server//CO2//Hete//ver//Plot//Processed_data//')
    radial_histogram('F://Server//CO2//Hete//ver//Plot//Histagram//')

def plot_data(file, title):
    #plt_file = pd.read_csv(Data_path + 'Processed_data//' + os.listdir(Data_path + 'Processed_data')[2])
    plt_file = file
    plt.figure(figsize=(12, 9))
    x = plt_file['Time']
    x = x*200
    y1 = plt_file['CO2 saturation']
    y2 = plt_file['Residual CO2 saturation']
    #y1_smooth = itp(x, y1, kind = 'nearest')(x)
    #y2_smooth = itp(x, y2, kind = 'nearest')(x)
    
    ax1 = plt.gca()
    ax1.tick_params(axis='both', labelsize=24)
    #ax1.scatter(x, y1)
    ax1.plot(x, y1, color='red', linewidth=5)
    
    for label in ax1.get_xticklabels() + ax1.get_yticklabels():
        label.set_fontweight('bold')
        label.set_fontname('Times New Roman')
        
    for label in ax1.get_yticklabels():
        label.set_color('red')
    
    ax2 = ax1.twinx()
    ax2.tick_params(axis='both', labelsize=24)
    #ax2.scatter(x, y2)
    ax2.plot(x, y2, color='blue', linewidth=5)
    
    for label in ax2.get_yticklabels():
        label.set_fontweight('bold')
        label.set_fontname('Times New Roman')
        label.set_color('blue')
    
    ax1.set_xlabel('Dimensionless time t*', fontdict={'name': 'Times New Roman', 'size': 28}, fontweight = 'bold', labelpad = 15)
    ax1.set_ylabel('Total CO2', fontdict={'name': 'Times New Roman', 'size': 28}, color = 'red', fontweight = 'bold', labelpad = 15)
    ax1.set_ylim(-0.018, 0.9)
    ax2.set_ylabel('Residual trapped CO2', fontdict={'name': 'Times New Roman', 'size': 28}, color = 'blue', fontweight = 'bold', labelpad = 15)
    ax2.set_ylim(-0.006, 0.3)
    ax1.legend()
    ax2.legend()
    plt.title(title)
    plt.show()
        
def get_xy(img, start_x=None, start_y=None, end_x=None, end_y=None):
    CO2 = img[:, :, 2]
    #start_x, start_y, end_x, end_y = 0, 0, len(CO2), len(CO2[0])
    for a in range(len(CO2)):
        if any (CO2[a,:]!=82):
            start_y = a
            break
        else:
            continue
    for b in range(len(CO2[0])):
        if any (CO2[:,b]!=82):
            start_x = b
            break
        else:
            continue
    for c in range(start_y, len(CO2)):
        if any (CO2[c,:]!=82):
            continue
        else:
            end_y = c
            break
    for d in range(start_x, len(CO2[0])):
        if any (CO2[:,d]!=82):
            continue
        else:
            end_x = d
            break
    return start_x, start_y, end_x, end_y

def img_prepro(path, sec_dir):
    imgs_CO2 = [file for file in os.listdir(path) if 'png' in file and 'CO2' in file]
    imgs_p = [file for file in os.listdir(path) if 'png' in file and not 'CO2' in file]
        
    for i in range(len(imgs_CO2)):
        im_CO2 = cv2.imread(path + imgs_CO2[i])
        im_p = cv2.imread(path + imgs_p[i])
        
        if sec_dir=='LHM' or sec_dir=='HLM' or sec_dir=='HML' in path:
            im_flip_CO2 = cv2.imread(path + imgs_CO2[i])
            im_flip_p = cv2.imread(path + imgs_p[i])
            im_flip_CO2 = cv2.flip(im_flip_CO2, 0)
            im_flip_p = cv2.flip(im_flip_p, 0)
            cv2.imwrite(path + imgs_CO2[i], im_flip_CO2)
            cv2.imwrite(path + imgs_p[i], im_flip_p)
        
        if len(im_CO2) >= 370 and len(im_CO2[0]) >= 130:
            start_x, start_y, end_x, end_y = get_xy(im_CO2)
            cropped_im_CO2 = im_CO2[start_y:end_y, start_x:end_x]
            cropped_im_p = im_p[start_y:end_y, start_x:end_x]
            cv2.imwrite(path + imgs_CO2[i], cropped_im_CO2)
            cv2.imwrite(path + imgs_p[i], cropped_im_p)
                        
    return imgs_CO2, imgs_p

#--------------------------------------------------------------------------------#

def pressure_ana(path, case_name,row,col,porosity,proc_file, mode='',area_num=45, time_split=2):
    
    csv_count = sum(1 for file in os.listdir(path) if file.endswith('.csv'))
    lists = [[] for _ in range(area_num)]
    p0_list = []

    for i in range(1, csv_count + 1):
        p = path + '//' + 'Data_{}.csv'.format(round(i*0.0001, 4) if i != 0 else 0)
        df_p = pd.read_csv(p)
        sorted_df = df_p.sort_values(by='X')
        if i % time_split == 0:
            p0 = sorted_df['p'].iloc[20]
            p0_list.append(p0)
                    
        seg_len = len(sorted_df) // area_num
        for k in range(area_num):
            start = k * seg_len
            end = (k+1) * seg_len if k < (area_num-1) else len(sorted_df)
            area_ave = sorted_df['p'][start:end].mean()
            if area_ave>1e7:
                area_ave = area_ave - 1e7
            lists[k].append(area_ave)
            
    '''time = [(t+1)*0.02 for t in range(100)]
    data_out = {'Time': time}
    for z in range(area_num):
        data_out['A{}'.format(z)] = lists[z]
    df = pd.DataFrame(data_out) 
    df.to_csv('F://Server//CO2//Hete//ver//Plot//Seg_pressure//' + case_name + '.csv', index=False)'''
    plt.figure(figsize = (12, 9))
    
    if mode=='time':
        time = [(t+1)*0.02 for t in range(100)]
        data_out = {'Time': time}
        for z in range(area_num):
            data_out['A{}'.format(z)] = lists[z]
        df = pd.DataFrame(data_out) 
        df.to_csv('F://Server//CO2//Hete//ver//Plot//Seg_pressure//' + case_name + '.csv', index=False)
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
        time = [(t+1)*0.02 for t in range(100)]
        data_out = {'Time': time}
        for z in range(area_num):
            data_out['A{}'.format(z)] = lists[z]
        df = pd.DataFrame(data_out) 
        df.to_csv('F://Server//CO2//Hete//ver//Plot//Seg_pressure//' + case_name + '.csv', index=False)
        x = []
        for t in range(area_num):
            x.append(t*(450/area_num) + (450/area_num)/2)
        
        sel_str = pd.read_csv('F://Server//CO2//Hete//ver//Plot//Curve//Img_selection.csv', header=None).iloc[row+1, col+1]
        sel_spli = sel_str.split(',')
        y0 = df.iloc[9, 1:]/1000
        y1 = df.iloc[eval(sel_spli[0])-1, 1:]/1000
        y2 = df.iloc[eval(sel_spli[1])-1, 1:]/1000
        y3 = df.iloc[eval(sel_spli[2])-1, 1:]/1000
        y4 = df.iloc[eval(sel_spli[3])-1, 1:]/1000  
        y5 = df.iloc[99, 1:]/1000
            
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
        ax.set_ylabel('ΔP (kPa)', fontdict={'name': 'Times New Roman', 'size': 28}, fontweight = 'bold', labelpad = 15)
        #ax.set_yscale('log')
        #ax.set_ylim(0,380)
        ax.xaxis.set_major_locator(MultipleLocator(50))
        ax.yaxis.set_major_locator(MultipleLocator(40))
        ax.legend(['t* = 0.2', 't* = {}'.format(round(eval(sel_spli[0])/50, 2)), 't* = {}'.format(round(eval(sel_spli[1])/50, 2)), 't* = {}'.format(round(eval(sel_spli[2])/50, 2)), 't* = {}'.format(round(eval(sel_spli[3])/50, 2)), 't* = 2'], prop={'family': 'Times New Roman', 'size': 20})
        plt.title(case_name)
        
    if mode=='energy':
        df_CO2_sat = pd.read_csv(proc_file)
        time = [(t*time_split)*0.02 for t in range(int(100/time_split))]
        data_out = {'Time': time}
        Ev = [z*100/(4.5*porosity*1e6) for z in p0_list]  #Mw/m3
        Em = []
        for s in range (len(p0_list)):
            if s<=50:
                Em.append(Ev[s]*1000/(720*df_CO2_sat['CO2 saturation'].iloc[s]))
            else:
                Em.append(Ev[s]*1000/(720*(df_CO2_sat['CO2 saturation'].iloc[s]-df_CO2_sat['CO2 saturation'].iloc[s-4:s].mean())))
        data_out['A1'] = Ev  #w/m^3
        data_out['A2'] = Em  #w/kgCO2
        df = pd.DataFrame(data_out) 
        df.to_csv('F://Server//CO2//Hete//ver//Plot//Energy_consum//' + case_name + '.csv', index=False)
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
        ax2.set_ylabel('Em (kw/kgCO₂)', fontdict={'name': 'Times New Roman', 'size': 28}, color = 'green', fontweight = 'bold', labelpad = 15)
        ax2.set_yscale('log')
        ax1.set_ylim(0,15)
        ax2.set_ylim(8,5000)
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
                y_max = df1['Residual CO2 saturation'].max()
                y_conser = df1['Residual CO2 saturation'].iloc[-3:-1].mean()
                y_total = df1['CO2 saturation'].iloc[-1]
                max_conser[0].append(y_max)
                max_conser[1].append(y_conser)
                max_conser[2].append(y_total)
                Heter.append(file1.split('_')[0] + '_' + file1.split('_')[-1].split('.')[0])
                
        for file2 in os.listdir('F://Server//CO2//Hete//ver//Plot//Energy_consum//'):
            if probe[t] in file2:
                df2 = pd.read_csv('F://Server//CO2//Hete//ver//Plot//Energy_consum//'+file2)
                y_em = df2['A1'].mean()
                y_ev = df2['A2'].mean()
                max_conser[3].append(y_em)
                max_conser[4].append(y_ev)
                
        data_out = {'Heterogeneity': Heter}
        data_out['Max'] = max_conser[0]
        data_out['Conserved'] = max_conser[1]
        data_out['Total'] = max_conser[2]
        data_out['A1'] = max_conser[3]
        data_out['A2'] = max_conser[4]
        out = pd.DataFrame(data_out)
        out.to_csv('F://Server//CO2//Hete//ver//Plot//Histagram//' + 'radial_histagram_{}.csv'.format(probe[t]), index=False)

def radial_histogram(csv_path):
    files = ['radial_histagram_CA100.csv', 'radial_histagram_CA125.csv', 'radial_histagram_CA150.csv']
    for csv_file in files:
        data = pd.read_csv(os.path.join(csv_path, csv_file))
        
        # Assuming 9 unique labels for 9 bins
        labels = data['Heterogeneity'].unique()
        num_labels = len(labels)
        angles = np.linspace(0, 2 * np.pi, num_labels, endpoint=False)
        
        # Setting up the plot
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 10))
        fig, ax1 = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 10))
        fig, ax2 = plt.subplots(figsize=(16, 10))
        
        # Calculating width of each bar to cover the entire bin
        width = 2 * pi / num_labels
        
        # Plot 'Max' data
        max_values = data['Max']
        ax.bar(angles, max_values, width=width, color='blue', alpha=0.6, edgecolor='black')
        
        # Plot 'Conserved' data on top
        conserved_values = data['Conserved']
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
        
        # Plot Ev data
        max_values = data['A1']
        ax2.bar(angles, max_values, width=width, color='orange', alpha=0.6, edgecolor='black')
        
        # Plot Em data on top
        conserved_values = data['A2']
        ax2.bar(angles, conserved_values, width=width, color='green', alpha=0.6, edgecolor='black')
        
        # Adjust grid and labels
        ax2.set_xticks(angles)  # Set ticks to the angles of the bars
        ax2.set_xticklabels(labels, fontdict={'name': 'Times New Roman', 'size': 20}) 
        #ax.set_xticklabels([])  # Remove tick labels
        ax2.set_yticklabels([])
        ax2.tick_params(axis='y', which='major', labelsize=24)
        ax2.yaxis.set_major_locator(MultipleLocator(20))
        ax2.set_yscale('log')
        ax2.set_ylim(0,100)
        
        # Optional: Customize grid and axis for better visibility and alignment
        ax2.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
        plt.title(csv_file)
        
        '''
        Ex = ['LMH_CA125','LMH_CA125g','LHM_CA125g','MHL_CA125','MHL_CA125g','MLH_CA125','MLH_CA125g','HML_CA125g','HLM_CA125g']
        Em_value = data['A1']
        Ev_value = data['A2']
        
        ind = np.arange(len(Em_value))
        ax2 = plt.gca()
        ax2.tick_params(axis='y', labelsize=24)
        #ax1.scatter(x, y1)
        ax2.bar(ind-0.2, Em_value, color='orange', width=0.4, alpha=0.85)
        
        for label in ax2.get_xticklabels() + ax2.get_yticklabels():
            label.set_fontweight('bold')
            label.set_fontname('Times New Roman')
            
        for label in ax2.get_yticklabels():
            label.set_color('orange')
        
        ax3 = ax2.twinx()
        ax3.tick_params(axis='y', labelsize=24)
        #ax2.scatter(x, y2)
        ax3.bar(ind+0.2, Ev_value, color='green', width=0.4, alpha=0.85)
        
        for label in ax3.get_yticklabels():
            label.set_fontweight('bold')
            label.set_fontname('Times New Roman')
            label.set_color('green')
        
        ax2.set_xticks(ind)
        ax2.set_xticklabels(Ex, fontdict={'name': 'Times New Roman', 'size': 18}, fontweight = 'bold', rotation=45)
        #ax3.set_xlabel(x, fontdict={'name': 'Times New Roman', 'size': 6}, fontweight = 'bold', labelpad = 15)
        ax2.set_ylabel('Ev (Mw/m³)', fontdict={'name': 'Times New Roman', 'size': 28}, color = 'orange', fontweight = 'bold', labelpad = 15)
        ax2.set_ylim(5, 10)
        #ax2.set_yscale('log')
        ax3.set_ylabel('Em (Mw/kgCO₂)', fontdict={'name': 'Times New Roman', 'size': 28}, color = 'green', fontweight = 'bold', labelpad = 15)
        #ax3.set_yscale('log')
        ax3.set_ylim(140, 175)
        #ax2.legend()
        #ax3.legend()
        plt.title(csv_file)'''
    plt.show()

#--------------------------------------------------------------------------------#

run_process()   
#get_max_conser('F://Server//CO2//Hete//ver//Plot//Processed_data//')
#radial_histogram('F://Server//CO2//Hete//ver//Plot//Histagram//')

'''case_name = ['hor_HML_CA150', 'hor_MHL_CA125', 'ver_LMH_CA150', 'ver_MLH_CA100']
path = ['E://Server//Hete//hor//HML//CO2_aquifer_immisible_CA150', 'E://Server//Hete//hor//MHL//CO2_aquifer_immisible_CA125', 'E://Server//Hete//ver//LMH//CO2_aquifer_immisible_CA150', 'E://Server//Hete//ver//MLH//CO2_aquifer_immisible_CA100']
for t in range(len(case_name)):
    pressure_ana(path[t], case_name[t], area_num=6)'''
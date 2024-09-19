#Setup required packages------------
import streamlit as st
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#import extra_streamlit_components as stx
import datetime
#Imports----------
from streamlit_option_menu import option_menu
from pandas import Series, DataFrame
from scipy.signal import find_peaks as fp
from scipy.signal import savgol_filter
from scipy.signal import peak_widths
from peakutils import indexes
from peakutils import baseline
from scipy import integrate
#----------------
#Title of the app
st.set_page_config(page_icon=None, layout="wide", initial_sidebar_state='auto')
st.title("4 Biomarker QIRT-ELISA Data Analysis")
#-----------------
#Sidebar
#with st.sidebar:
#    selected = option_menu("Steps", ["Data Loading", 'Signal Smoothing', 'Baseline Correction', 'Peak Analysis'], 
#        icons=['1-circle-fill', '2-circle-fill', '3-circle-fill', '4-circle-fill'], menu_icon="list", default_index=0)
#----------------
#st.sidebar.markdown('''
# Sections
#- [Data Loading](#section-1)
#- [Signal Smoothing](#section-2)
#- [Baseline Correction](#section-3)
#- [Peak Analysis](#section-4)                                        
#''', unsafe_allow_html=True)    
#----------------------
#Setup tabs
tab1, tab2, tab3, tab4, tab5 = st.tabs(["Data Loading", "Signal Smoothing", "Baseline Correction", "Peak Analysis", "Results"])
#----------------------
#Load data
with tab1:
    st.header('Load Your Data')
    #@st.cache_data
    #st.cache_data.clear()
    def load_data(file):
     data = pd.read_csv(file, usecols=['#Digilent WaveForms Oscilloscope Acquisition','Unnamed: 1','Unnamed: 2', 'Unnamed: 3', 'Unnamed: 4']).rename(columns = {'Unnamed: 1': 'Signal605', 'Unnamed: 4': 'Signal655', 'Unnamed: 3': 'Signal525', 'Unnamed: 2': 'Signal565',  '#Digilent WaveForms Oscilloscope Acquisition': 'time'})
     return data

    uploaded_file = st.file_uploader ("upload your data file", type={"csv", "txt"})

    if uploaded_file is None:
        st.info(" Please upload a file", icon="ℹ️")
        st.stop()

    df = load_data(uploaded_file)
 #------------------
 #data clean-up to remove the headers from the csv files
    df.drop(df.index[pd.Series(range(0,15))],inplace=True)
    df=df.astype(float)
    #df['Signal655'] = df['Signal655'] * 1
    with st.expander("Data Preview"):
            st.dataframe(df)
 #------------------
 #Save the name of your timepoint
 #    experiment_replicate = st.text_input(
 #    "Enter the sample timepoint - replicate:", placeholder="Formatted as t##-r##"
 #    )
 #    experiment_date = st.date_input("Enter the date of experiment", value=None)
 #------------------
 #Select the QDot channels you want to analyze
 #    channels = st.multiselect(
 #    'Which QDot channel do you want to analyze?', 
 #    ['QDot605', 'QDot655'], help = 'Please select at least one Qdot channel')
 #    if channels is None:
 #        st.info(" Please choose your channel/s", icon="ℹ️")
 #        st.stop()

 #    for item in channels:
 #        if item == 'QDot605':            
            
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        multiplier1 = st.number_input('Preamplifier sensitivity for QDot525:', min_value=1, help = 'Please enter in pA')
        df['Signal525'] = df['Signal525'] * multiplier1 
    with col2:
        multiplier2 = st.number_input('Preamplifier sensitivity for QDot565:', min_value=1, help = 'Please enter in pA')
        df['Signal565'] = df['Signal565'] * multiplier2
    with col3:
        multiplier3 = st.number_input('Preamplifier sensitivity for QDot605:', min_value=1, help = 'Please enter in pA')
        df['Signal605'] = df['Signal605'] * multiplier3
    with col4:
        multiplier4 = st.number_input('Preamplifier sensitivity for QDot655:', min_value=1, help = 'Please enter in pA')
        df['Signal655'] = df['Signal655'] * multiplier4    

    if multiplier1 is 1:
        st.info(" Please enter the sensitivities values", icon="ℹ️")
        st.stop()
    if multiplier2 is 1:
        st.info(" Please enter the sensitivities values", icon="ℹ️")
        st.stop()
    if multiplier3 is 1:
        st.info(" Please enter the sensitivities values", icon="ℹ️")
        st.stop() 
    if multiplier4 is 1:
        st.info(" Please enter the sensitivities values", icon="ℹ️")
        st.stop()      


    st.sidebar.info('Move to the *Signal Smoothing* tab once you have entered your preamplifier sensitivities.', icon="1️⃣")
    with st.expander("Data Preview"):
            st.dataframe(df)
with tab2:
    st.header('Signal Smoothing')
    st.write('By default, there is no signal smoothing. Please check the boxes below to apply a Savitzky-Golay smoothing filter if necessary.',icon="1️⃣")
    #st.info('The parameters to choose for the filter is: *the length of the filter window* and the order of the polynomial used to fit the data')
    
    col1, col2, col3, col4 = st.columns(4) 
    with col1:
        st.header("QDot 525")
    
        filter_window1 = 1
        order_value1 = 0
        Signal525_smoothed = savgol_filter(df.Signal525, filter_window1, order_value1, deriv=0, delta=1.0, axis=-1, mode='interp', cval=0.0)
        
        fig_smooth_525 = plt.figure(figsize=(9, 7))
        sns.lineplot(y =Signal525_smoothed, x = df.time, label = "QDot 525 Signal")
        plt.xlabel("Time [sec]")
        plt.ylabel("Fluorscent Signal [pA]")
        plt.title("Signals")
        st.pyplot(fig_smooth_525)

        change_smooth_525 = st.checkbox("Apply a signal smoothing filter for QDot 525")
        if change_smooth_525:
            filter_window1 = st.number_input("Length of the filter window", min_value=1  , max_value=None, value=None, placeholder="Type a number...", key='fw1', help = 'Please enter a value less the length of your data frame. For no signal smoothing, enter 1')
            if filter_window1 is None:
                st.info(" Please choose a filter window", icon="ℹ️")
                st.stop() 
            order_value1 = st.number_input("Order of the polynomial", min_value=0, max_value=filter_window1, value=None, placeholder="Type a number...", key='ov1', help = 'Please enter a value less the filter window. For no signal smoothing, enter 0')
            if order_value1 is None:
                st.info(" Please choose the order of the filter", icon="ℹ️")
                st.stop() 
            Signal525_smoothed = savgol_filter(df.Signal525, filter_window1, order_value1, deriv=0, delta=1.0, axis=-1, mode='interp', cval=0.0)
        
            fig_smooth_525 = plt.figure(figsize=(9, 7))
            sns.lineplot(y =df.Signal525, x = df.time, label = "525 Raw Signal") 
            sns.lineplot(y =Signal525_smoothed, x = df.time, label = "QDot 525 Signal")
            plt.xlabel("Time [sec]")
            plt.ylabel("Fluorscent Signal [pA]")
            plt.title("Signals")
            st.pyplot(fig_smooth_525)

    with col2:
        st.header("QDot 565")
            
        filter_window2 = 1
        order_value2 = 0 
        Signal565_smoothed = savgol_filter(df.Signal565, filter_window2, order_value2, deriv=0, delta=1.0, axis=-1, mode='interp', cval=0.0)

        fig_smooth_565 = plt.figure(figsize=(9, 7))   
        sns.lineplot(y =Signal565_smoothed, x = df.time, label = "QDot 565 Signal")
        plt.title("Signals")
        plt.xlabel("Time [sec]")
        plt.ylabel("Fluorscent Signal [pA]")
        st.pyplot(fig_smooth_565)
    
        st.sidebar.info('Move to the *Baseline correction* tab if you are satisfied with your filter selection.', icon="2️⃣")

        change_smooth_565 = st.checkbox("Apply a signal smoothing filter for QDot 565")
        if change_smooth_565:
            filter_window2 = st.number_input("Length of the filter window", min_value=1  , max_value=None, value=None, placeholder="Type a number...", key='fw2', help = 'Please enter a value less the length of your data frame. For no signal smoothing, enter 1')
            if filter_window2 is None:
                st.info(" Please choose a filter window", icon="ℹ️")
                st.stop() 
            order_value2 = st.number_input("Order of the polynomial", min_value=0, max_value=filter_window2, value=None, placeholder="Type a number...", key='ov2', help = 'Please enter a value less the filter window. For no signal smoothing, enter 0')
            if order_value2 is None:
                st.info(" Please choose the order of the filter", icon="ℹ️")
                st.stop()  
            Signal565_smoothed = savgol_filter(df.Signal565, filter_window2, order_value2, deriv=0, delta=1.0, axis=-1, mode='interp', cval=0.0)
            fig_smooth_565 = plt.figure(figsize=(9, 7))   
            sns.lineplot(y =df.Signal565, x = df.time, label = "565 Raw Signal")
            sns.lineplot(y =Signal565_smoothed, x = df.time, label = "QDot 565 Signal")
            plt.title("Signals")
            plt.xlabel("Time [sec]")
            plt.ylabel("Fluorscent Signal [pA]")
            st.pyplot(fig_smooth_565) 
    
    with col3:
        st.header("QDot 605")
    
        filter_window3 = 1
        order_value3 = 0
        Signal605_smoothed = savgol_filter(df.Signal605, filter_window3, order_value3, deriv=0, delta=1.0, axis=-1, mode='interp', cval=0.0)

        fig_smooth_605 = plt.figure(figsize=(9, 7))
        sns.lineplot(y =Signal605_smoothed, x = df.time, label = "QDot 605 Signal")
        plt.xlabel("Time [sec]")
        plt.ylabel("Fluorscent Signal [pA]")
        plt.title("Signals")
        st.pyplot(fig_smooth_605)

        change_smooth_605 = st.checkbox("Apply a signal smoothing filter for QDot 605")
        if change_smooth_605:
            filter_window3 = st.number_input("Length of the filter window", min_value=1  , max_value=None, value=None, placeholder="Type a number...", key='fw3', help = 'Please enter a value less the length of your data frame. For no signal smoothing, enter 1')
            if filter_window3 is None:
                st.info(" Please choose a filter window", icon="ℹ️")
                st.stop() 
            order_value3 = st.number_input("Order of the polynomial", min_value=0, max_value=filter_window3, value=None, placeholder="Type a number...", key='ov3', help = 'Please enter a value less the filter window. For no signal smoothing, enter 0')
            if order_value3 is None:
                st.info(" Please choose the order of the filter", icon="ℹ️")
                st.stop() 
            Signal605_smoothed = savgol_filter(df.Signal605, filter_window3, order_value3, deriv=0, delta=1.0, axis=-1, mode='interp', cval=0.0)
        
            fig_smooth_605 = plt.figure(figsize=(9, 7))
            sns.lineplot(y =df.Signal605, x = df.time, label = "605 Raw Signal") 
            sns.lineplot(y =Signal605_smoothed, x = df.time, label = "QDot 605 Signal")
            plt.xlabel("Time [sec]")
            plt.ylabel("Fluorscent Signal [pA]")
            plt.title("Signals")
            st.pyplot(fig_smooth_605)

    with col4:
        st.header("QDot 655")
    
        filter_window4 = 1
        order_value4 = 0
        Signal655_smoothed = savgol_filter(df.Signal655, filter_window4, order_value4, deriv=0, delta=1.0, axis=-1, mode='interp', cval=0.0)
        
        fig_smooth_655 = plt.figure(figsize=(9, 7))
        sns.lineplot(y =Signal655_smoothed, x = df.time, label = "QDot 655 Signal")
        plt.xlabel("Time [sec]")
        plt.ylabel("Fluorscent Signal [pA]")
        plt.title("Signals")
        st.pyplot(fig_smooth_655)

        change_smooth_655 = st.checkbox("Apply a signal smoothing filter for QDot 655")
        if change_smooth_655:
            filter_window4 = st.number_input("Length of the filter window", min_value=1  , max_value=None, value=None, placeholder="Type a number...", key='fw4', help = 'Please enter a value less the length of your data frame. For no signal smoothing, enter 1')
            if filter_window4 is None:
                st.info(" Please choose a filter window", icon="ℹ️")
                st.stop() 
            order_value4 = st.number_input("Order of the polynomial", min_value=0, max_value=filter_window4, value=None, placeholder="Type a number...", key='ov4', help = 'Please enter a value less the filter window. For no signal smoothing, enter 0')
            if order_value4 is None:
                st.info(" Please choose the order of the filter", icon="ℹ️")
                st.stop() 
            Signal655_smoothed = savgol_filter(df.Signal655, filter_window4, order_value4, deriv=0, delta=1.0, axis=-1, mode='interp', cval=0.0)
        
            fig_smooth_655 = plt.figure(figsize=(9, 7))
            sns.lineplot(y =df.Signal655, x = df.time, label = "655 Raw Signal") 
            sns.lineplot(y =Signal655_smoothed, x = df.time, label = "QDot 655 Signal")
            plt.xlabel("Time [sec]")
            plt.ylabel("Fluorscent Signal [pA]")
            plt.title("Signals")
            st.pyplot(fig_smooth_655)


with tab3:
    
    st.header('Baseline Correction')
    st.write('Remove the background noise')

    col1, col2, col3, col4 = st.columns(4)

    with col1:
     
     base_deg1 = 100

     base_line1 = baseline(Signal525_smoothed, deg = base_deg1) #Create the baseline
     Signal525_smoothed_adjusted = Signal525_smoothed - base_line1 #Subtract the baseline

     fig_baseline_525 = plt.figure(figsize=(9, 7))
     sns.lineplot(y = base_line1, x = df.time, label = "Generated Baseline for QDot 525")
     sns.lineplot(y = Signal525_smoothed, x = df.time, label = "Adjusted QDot 605 Signal")
     plt.title( "Signals")
     plt.xlabel("Time [sec]")
     plt.ylabel("Fluorscent Signal [pA]")
     st.pyplot(fig_baseline_525)

    with col2:

     base_deg2 = 100

     base_line2 = baseline(Signal565_smoothed, deg = base_deg2) #Create the baseline
     Signal565_smoothed_adjusted = Signal565_smoothed - base_line2 #Subtract the baseline

     fig_baseline_565 = plt.figure(figsize=(9, 7))
     sns.lineplot(y = base_line2, x = df.time, label = "Generated Baseline for QDot 565")
     sns.lineplot(y = Signal565_smoothed, x = df.time, label = "Adjusted QDot 565 Signal")
     plt.title( "Signals")
     plt.xlabel("Time [sec]")
     plt.ylabel("Fluorscent Signal [pA]")
     st.pyplot(fig_baseline_565)

    with col3:
     
     base_deg3 = 100

     base_line3 = baseline(Signal605_smoothed, deg = base_deg1) #Create the baseline
     Signal605_smoothed_adjusted = Signal605_smoothed - base_line3 #Subtract the baseline

     fig_baseline_605 = plt.figure(figsize=(9, 7))
     sns.lineplot(y = base_line3, x = df.time, label = "Generated Baseline for QDot 605")
     sns.lineplot(y = Signal605_smoothed, x = df.time, label = "Adjusted QDot 605 Signal")
     plt.title( "Signals")
     plt.xlabel("Time [sec]")
     plt.ylabel("Fluorscent Signal [pA]")
     st.pyplot(fig_baseline_605)

    with col4:
     
     base_deg4 = 100

     #base_deg1 = st.number_input("Degree of the baseline for QDot 605", min_value=1  , max_value=None, value=None, placeholder="Type a number...")
     base_line4 = baseline(Signal655_smoothed, deg = base_deg1) #Create the baseline
     Signal655_smoothed_adjusted = Signal655_smoothed - base_line4 #Subtract the baseline

     fig_baseline_655 = plt.figure(figsize=(9, 7))
     sns.lineplot(y = base_line4, x = df.time, label = "Generated Baseline for QDot 655")
     sns.lineplot(y = Signal655_smoothed, x = df.time, label = "Adjusted QDot 655 Signal")
     plt.title( "Signals")
     plt.xlabel("Time [sec]")
     plt.ylabel("Fluorscent Signal [pA]")
     st.pyplot(fig_baseline_655)

     st.sidebar.info('Move to the *Peak Analysis* tab after the signals are adjusted.', icon="3️⃣")
    

with tab4:
    st.header('Peak Analysis')
    #build new data frames for peak finding
    df21 = pd.DataFrame({'time': df.time,'Signal':Signal525_smoothed_adjusted}).set_index('time')
    df22 = pd.DataFrame({'time': df.time,'Signal':Signal565_smoothed_adjusted}).set_index('time')
    df23 = pd.DataFrame({'time': df.time,'Signal':Signal605_smoothed_adjusted}).set_index('time')
    df24 = pd.DataFrame({'time': df.time,'Signal':Signal655_smoothed_adjusted}).set_index('time') 
    st.write('Find the peaks in the adjusted signals and integrate the found peak.')


    #def set_inputs_zero():
    #    st.session_state.h2 = 0
    #    st.session_state.prom2 = 0
    
        
    col1, col2, col3, col4 = st.columns(4)

    with col1:
         st.header("QDot 525")
         
         maxval_525 = df21.Signal.max()
         h1 = 0.1*maxval_525
         prom1 = h1/2

         def set_inputs_zero():
          st.session_state.h2 = 0
          st.session_state.prom2 = 0
  
         peaks21, property1 = fp(df21.Signal,
            height = h1,
            prominence = prom1,
            distance = 10, wlen = 350)
         peak_found_525 = property1


         change525 = st.checkbox("Change the pre-set peaks for QDot 525")
         if change525:
            new_prom1 = st.slider(label = "Choose the new prominence", min_value = 1 , max_value = None, value = None, step = 1, key = "changeprom1")
            prom1 = new_prom1
            peaks21, property1 = fp(df21.Signal,
               height = h1,
               prominence = prom1,
               distance = 1, wlen = 350)
            peak_found_605 = property1
         
         df21.plot()
         sns.scatterplot(data = df21.iloc[peaks21], x = 'time', y = 'Signal', 
                color = 'red', alpha = 0.5)
         st.pyplot(plt.gcf( ))
            
    with col2:
         st.header("QDot 565")

         maxval_565 = df22.Signal.max()
         h2 = 0.1*maxval_565
         prom2 = h2/2

   
         peaks22, property2 = fp(df22.Signal,
            height = h2,
            prominence = prom2,
            distance = 10, wlen = 350)
         peak_found_565 = property2



         change565 = st.checkbox("Change the preset peaks for QDot 565")
         if change565:
            new_prom2 = st.slider(label = "Choose the new prominence", min_value = 1 , max_value = None, value = None, step = 1, key = "changeprom2")
            prom2 = new_prom2
            peaks22, property2 = fp(df22.Signal,
               height = h2,
               prominence = prom2,
               distance = 10, wlen = 350)
            peak_found_565 = property2


         df22.plot()
         sns.scatterplot(data = df22.iloc[peaks22], x = 'time', y = 'Signal', 
                color = 'red', alpha = 0.5)
         st.pyplot(plt.gcf( ))
         st.sidebar.info('Move to the *Results* tab if you are satisfied with the identified peaks.', icon="4️⃣")
    
    with col3:
         st.header("QDot 605")

         maxval_605 = df23.Signal.max()
         h3 = 0.1*maxval_605
         prom3 = h3/2

   
         peaks23, property3 = fp(df23.Signal,
            height = h3,
            prominence = prom3,
            distance = 10, wlen = 350)
         peak_found_605 = property3



         change605 = st.checkbox("Change the preset peaks for QDot 605")
         if change605:
            new_prom3 = st.slider(label = "Choose the new prominence", min_value = 1 , max_value = None, value = None, step = 1, key = "changeprom3")
            prom3 = new_prom3
            peaks23, property3 = fp(df23.Signal,
               height = h3,
               prominence = prom3,
               distance = 10, wlen = 350)
            peak_found_605 = property3


         df23.plot()
         sns.scatterplot(data = df23.iloc[peaks23], x = 'time', y = 'Signal', 
                color = 'red', alpha = 0.5)
         st.pyplot(plt.gcf( ))
         st.sidebar.info('Move to the *Results* tab if you are satisfied with the identified peaks.', icon="4️⃣")
    
    with col4:
         st.header("QDot 655")

         maxval_655 = df24.Signal.max()
         h4 = 0.1*maxval_655
         prom4 = h4/2

   
         peaks24, property4 = fp(df24.Signal,
            height = h4,
            prominence = prom4,
            distance = 10, wlen = 350)
         peak_found_655 = property4



         change655 = st.checkbox("Change the preset prominence for QDot 655")
         if change655:
            new_prom4 = st.slider(label = "Choose the new prominence", min_value = 1 , max_value = None, value = None, step = 1, key = "changeprom4")
            prom4 = new_prom4
            peaks24, property4 = fp(df24.Signal,
               height = h4,
               prominence = prom4,
               distance = 10, wlen = 350)
            peak_found_655 = property4


         df24.plot()
         sns.scatterplot(data = df24.iloc[peaks24], x = 'time', y = 'Signal', 
                color = 'red', alpha = 0.5)
         st.pyplot(plt.gcf( ))
         st.sidebar.info('Move to the *Results* tab if you are satisfied with the identified peaks.', icon="4️⃣")
                       
         

with tab5:
    st.header('Results')
        
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        df21_2 = pd.DataFrame({'time': df.time,'Signal':Signal525_smoothed_adjusted}).set_index('time')
        df21_3 = pd.DataFrame({'time': df.time,'Signal':Signal525_smoothed_adjusted})
        time_array_525 = df21_3[["time"]].to_numpy( )
        peak_number_525 = len(peak_found_525['peak_heights'])
        n_525 = 0
        area_idv_525 = 0
        area_sum_525 = 0
        area_peaks_525 = 0 
        array_525_peaks = np.array([])

        for n_525 in range(peak_number_525):
         left_base_index_525 = peak_found_525['left_bases'][n_525]
         right_base_index_525 = peak_found_525['right_bases'][n_525]
         left_base_525 = time_array_525[left_base_index_525][0]
         right_base_525 = time_array_525[right_base_index_525][0]

         df5_525 = df21_2.Signal[left_base_525:right_base_525]
         timestamp2_525 = df5_525.index
         area_idv_525 = np.trapz(df21_2.Signal[left_base_525:right_base_525], x = timestamp2_525)
         array_525_peaks = np.append(array_525_peaks, area_idv_525)
         area_sum_525 = area_sum_525 + area_idv_525 
         n_525 = n_525 + 1
         area_peaks_525 = area_sum_525/peak_number_525

        y525 = [1.14493664, -0.5470587, -0.597878, 10.2377709, 11.0072546, 7.76470096, 22.4991946, 24.8191155, 24.1932251, 40.8335942, 41.9158234, 38.4410576, 53.0854621, 61.1917111, 70.4062313]
        x525 = [0, 0, 0, 10, 10, 10, 40, 40, 40, 80, 80, 80, 100, 100, 100]
        x_min_525 = min(x525)
        x_max_525 = max(x525)
        step525 = 0.01
        x_new_525 = np.arange(x_min_525, x_max_525, step525)
        #A1 = 0.2465
        B1_525 = 1.361
        #C1 = 27397571
        #D1 = 4208925
        m_525 = 0.5594
        #def log4pl_conc_gluc(y_exp_gluc):
        #    return C1*((A1-D1)/(y_exp_gluc-D1)-1)**(1/B1)
        #yfit1_new = ((A1-D1)/(1.0+((x_new_1/C1)**B1))) + D1

        def log4pl_conc_gluc(y_exp_glucose):
            return (y_exp_glucose-B1_525)/m_525 
        yfit1_new_525 = (m_525*x_new_525)+B1_525
        
        f_525 = area_peaks_525
        f_0_glucose= 0.812388566
        normalized_f_525 = (f_525 - f_0_glucose)/f_0_glucose
        y_exp_glucose = normalized_f_525


        glucose_concent = log4pl_conc_gluc(y_exp_glucose)
        fig_calib_glucose = plt.figure(figsize=(9, 7))
        plt.plot(x525, y525, 'r+', label="Experimental replicates")
        plt.plot(x_new_525, yfit1_new_525, label="Linear fit")
        plt.plot(glucose_concent,y_exp_glucose,'bd')
        plt.xlabel('Glucose Concentration [mM]')
        plt.ylabel('Fluorscent Signal (peak AUC)')
        plt.legend(loc='best', fancybox=True, shadow=True)
        st.pyplot(fig_calib_glucose)

        st.write('The individual peak AUCs are', array_525_peaks)
        st.write(df21.iloc[peaks21])
        st.write('The average QDOT525 peak AUC is:', area_peaks_525)
        st.write('The normalized average peak AUC is:', normalized_f_525)
        st.write('The corresponding glucose concentrations in [pM] is:', glucose_concent)

    with col2:
        df22_2 = pd.DataFrame({'time': df.time,'Signal':Signal565_smoothed_adjusted}).set_index('time')
        df22_3 = pd.DataFrame({'time': df.time,'Signal':Signal565_smoothed_adjusted})
        time_array_565 = df22_3[["time"]].to_numpy( )
        peak_number_565 = len(peak_found_565['peak_heights'])
        n_565 = 0
        area_idv_565 = 0
        area_sum_565 = 0
        area_peaks_565 = 0 
        array_565_peaks = np.array([])

        for n_565 in range(peak_number_565):
         left_base_index_565 = peak_found_565['left_bases'][n_565]
         right_base_index_565 = peak_found_565['right_bases'][n_565]
         left_base_565 = time_array_565[left_base_index_565][0]
         right_base_565 = time_array_565[right_base_index_565][0]

         df5_565 = df22_2.Signal[left_base_565:right_base_565]
         timestamp2_565 = df5_565.index
         area_idv_565 = np.trapz(df22_2.Signal[left_base_565:right_base_565], x = timestamp2_565)
         array_565_peaks = np.append(array_565_peaks, area_idv_565)
         area_sum_565 = area_sum_565 + area_idv_565 
         n_565 = n_565 + 1
         area_peaks_565 = area_sum_565/peak_number_565

        y565 = [1.14493664, -0.5470587, -0.597878, 10.2377709, 11.0072546, 7.76470096, 22.4991946, 24.8191155, 24.1932251, 40.8335942, 41.9158234, 38.4410576, 53.0854621, 61.1917111, 70.4062313]
        x565 = [0, 0, 0, 10, 10, 10, 40, 40, 40, 80, 80, 80, 100, 100, 100]
        x_min_565 = min(x565)
        x_max_565 = max(x565)
        step_565 = 0.01
        x_new_565 = np.arange(x_min_565, x_max_565, step_565)
        #A1 = 0.2465
        B1_565 = 1.361
        #C1 = 27397571
        #D1 = 4208925
        m_565 = 0.5594
        #def log4pl_conc_gluc(y_exp_gluc):
        #    return C1*((A1-D1)/(y_exp_gluc-D1)-1)**(1/B1)
        #yfit1_new = ((A1-D1)/(1.0+((x_new_1/C1)**B1))) + D1

        def log4pl_conc_cpep(y_exp_cpep):
            return (y_exp_cpep-B1_565)/m_565 
        yfit1_new_565 = (m_565*x_new_565)+B1_565
        
        f_565 = area_peaks_565
        f_0_cpep = 0.812388566
        normalized_f_565 = (f_565 - f_0_cpep)/f_0_cpep
        y_exp_cpep = normalized_f_565


        cpep_concent = log4pl_conc_cpep(y_exp_cpep)
        fig_calib_cpep = plt.figure(figsize=(9, 7))
        plt.plot(x565, y565, 'r+', label="Experimental replicates")
        plt.plot(x_new_565, yfit1_new_565, label="Linear fit")
        plt.plot(cpep_concent,y_exp_cpep,'bd')
        plt.xlabel('C-peptide Concentration [pM]')
        plt.ylabel('Fluorscent Signal (peak AUC)')
        plt.legend(loc='best', fancybox=True, shadow=True)
        st.pyplot(fig_calib_cpep)

        st.write('The individual peak AUCs are', array_565_peaks)
        st.write(df22.iloc[peaks22])
        st.write('The average QDOT565 peak AUC is:', area_peaks_565)
        st.write('The normalized average peak AUC is:', normalized_f_565)
        st.write('The corresponding C-peptide concentrations in [pM] is:', cpep_concent)

    with col3:    


        df23_2 = pd.DataFrame({'time': df.time,'Signal':Signal605_smoothed_adjusted}).set_index('time')
        df23_3 = pd.DataFrame({'time': df.time,'Signal':Signal605_smoothed_adjusted})
        time_array_605 = df23_3[["time"]].to_numpy( )
        peak_number_605 = len(peak_found_605['peak_heights'])
        n_605 = 0
        area_idv_605 = 0
        area_sum_605 = 0
        area_peaks_605 = 0
        array_605_peaks = np.array([])

        for n_605 in range(peak_number_605):
           left_base_index_605 = peak_found_605['left_bases'][n_605]
           right_base_index_605 = peak_found_605['right_bases'][n_605]
           left_base_605 = time_array_605[left_base_index_605][0]
           right_base_605 = time_array_605[right_base_index_605][0]
           df5_605 = df23_2.Signal[left_base_605:right_base_605]
           timestamp2_605 = df5_605.index
           area_idv_605 = np.trapz(df23_2.Signal[left_base_605:right_base_605], x = timestamp2_605)
           array_605_peaks = np.append(array_605_peaks, area_idv_605)
           area_sum_605 = area_sum_605 + area_idv_605 
           n_605 = n_605 + 1
           area_peaks_605 = area_sum_605/peak_number_605


        y605 = [0.58123804, -0.281096, -0.3001421, 8.06791665,	9.20469324,	7.49843809, 17.1460606, 14.6611185, 15.8027406, 31.3812569, 38.9815661, 35.3228082, 57.5532343, 63.0998399, 62.60464]
        x605 = [0, 0, 0, 100, 100, 100, 200, 200, 200, 400, 400, 400, 1000, 1000, 1000]
        x_min_605 = min(x605)
        x_max_605 = max(x605)
        step_605 = 0.01
        x_new_605 = np.arange(x_min_605, x_max_605, step_605)
        A2_605 = 0.3460
        B2_605 = 1.495
        C2_605 = 500.5
        D2_605 = 82.79
        def log4pl_conc_glucagon(y_exp_glucagon):
            return C2_605*((A2_605-D2_605)/(y_exp_glucagon-D2_605)-1)**(1/B2_605)
        yfit2_605 = ((A2_605-D2_605)/(1.0+((x_new_605/C2_605)**B2_605))) + D2_605


        f_605 = area_peaks_605
        f_0_glucagon = 1.426343498
        normalized_f_605 = (f_605 - f_0_glucagon)/f_0_glucagon
        y_exp_glucagon = normalized_f_605


        #y_exp_insulin = area_peaks_655
        glucagon_concent = log4pl_conc_glucagon(y_exp_glucagon)
        fig_calib_glucagon = plt.figure(figsize=(9, 7))
        plt.plot(x605, y605, 'r+', label="Experimental replicates")
        plt.plot(x_new_605, yfit2_605, label="Sigmoidal, 4PL non-linear fit")
        plt.plot(glucagon_concent,y_exp_glucagon,'bd')
        plt.xlabel('Glucagon Concentration [pM]')
        plt.ylabel('Fluorscent Signal (peak AUC)')
        plt.legend(loc='best', fancybox=True, shadow=True)
        st.pyplot(fig_calib_glucagon)

        st.write('The individual peak AUCs are', array_605_peaks)
        st.write(df23.iloc[peaks23])
        st.write('The average QDot605 peak AUC is:', area_peaks_605)
        st.write('The normalized average peak AUC is:', normalized_f_605)
        st.write('The corresponding Glucagon concentrations in [pM] is:', glucagon_concent)

    with col4:
        df24_2 = pd.DataFrame({'time': df.time,'Signal':Signal655_smoothed_adjusted}).set_index('time')
        df24_3 = pd.DataFrame({'time': df.time,'Signal':Signal655_smoothed_adjusted})
        time_array_655 = df24_3[["time"]].to_numpy( )
        peak_number_655 = len(peak_found_655['peak_heights'])
        n_655 = 0
        area_idv_655 = 0
        area_sum_655 = 0
        area_peaks_655 = 0 
        array_655_peaks = np.array([])

        for n_655 in range(peak_number_655):
         left_base_index_655 = peak_found_655['left_bases'][n_655]
         right_base_index_655 = peak_found_655['right_bases'][n_655]
         left_base_655 = time_array_655[left_base_index_655][0]
         right_base_655 = time_array_655[right_base_index_655][0]

         df5_655 = df24_2.Signal[left_base_655:right_base_655]
         timestamp2_655 = df5_655.index
         area_idv_655 = np.trapz(df24_2.Signal[left_base_655:right_base_655], x = timestamp2_655)
         array_655_peaks = np.append(array_655_peaks, area_idv_655)
         area_sum_655 = area_sum_655 + area_idv_655 
         n_655 = n_655 + 1
         area_peaks_655 = area_sum_655/peak_number_655

        y655 = [1.14493664, -0.5470587, -0.597878, 10.2377709, 11.0072546, 7.76470096, 22.4991946, 24.8191155, 24.1932251, 40.8335942, 41.9158234, 38.4410576, 53.0854621, 61.1917111, 70.4062313]
        x655 = [0, 0, 0, 10, 10, 10, 40, 40, 40, 80, 80, 80, 100, 100, 100]
        x_min_655 = min(x655)
        x_max_655 = max(x655)
        step1_655 = 0.01
        x_new_655 = np.arange(x_min_655, x_max_655, step1_655)
        #A1 = 0.2465
        B1_655 = 1.361
        #C1 = 27397571
        #D1 = 4208925
        m_655 = 0.5594
        #def log4pl_conc_gluc(y_exp_gluc):
        #    return C1*((A1-D1)/(y_exp_gluc-D1)-1)**(1/B1)
        #yfit1_new = ((A1-D1)/(1.0+((x_new_1/C1)**B1))) + D1

        def log4pl_conc_ins(y_exp_ins):
            return (y_exp_ins-B1_655)/m_655 
        yfit1_655 = (m_655*x_new_655)+B1_655
        
        f_655 = area_peaks_655
        f_0_ins = 0.812388566
        normalized_f_655 = (f_655 - f_0_ins)/f_0_ins
        y_exp_ins = normalized_f_655


        ins_concent = log4pl_conc_ins(y_exp_ins)
        fig_calib_ins = plt.figure(figsize=(9, 7))
        plt.plot(x655, y655, 'r+', label="Experimental replicates")
        plt.plot(x_new_655, yfit1_655, label="Linear fit")
        plt.plot(ins_concent,y_exp_ins,'bd')
        plt.xlabel('Insulin Concentration [pM]')
        plt.ylabel('Fluorscent Signal (peak AUC)')
        plt.legend(loc='best', fancybox=True, shadow=True)
        st.pyplot(fig_calib_ins)

        st.write('The individual peak AUCs are', array_655_peaks)
        st.write(df24.iloc[peaks24])
        st.write('The average QDot655 peak AUC is:', area_peaks_655)
        st.write('The normalized average peak AUC is:', normalized_f_655)
        st.write('The corresponding Insulin concentrations in [pM] is:', ins_concent)
       

        

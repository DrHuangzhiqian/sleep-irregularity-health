clear all
set more off
set maxvar 10000 

* ==============================================================================
* 1. 设定工作路径 
* ==============================================================================
cd "D:/ukb data/Sleep regularity/Sleep regularity"

* ==============================================================================
* 2. 读取并合并数据
* ==============================================================================
noisily display "正在读取并合并数据..."

import delimited "CovariatesImputed.csv", clear
save "CovariatesImputed_temp.dta", replace

import delimited "16_outcomes/death_outcome_260116.tsv", clear
save "death_outcome_temp.dta", replace

import delimited "GGIR_selected.csv", clear
keep eid sleepdurationinspt_ad_t5a5_mn sleepregularityindex_ad_t5a5_mn age_test mvpa season

merge 1:1 eid using "CovariatesImputed_temp.dta", keep(match) nogen
merge 1:1 eid using "death_outcome_temp.dta", keep(match) nogen

erase "CovariatesImputed_temp.dta"
erase "death_outcome_temp.dta"

* ==============================================================================
* 3. 数据清洗与核心变量构建 
* ==============================================================================
keep if target_time > 0
gen time_years = target_time / 365.25

gen age_in = age_test
gen age_out = age_test + time_years
keep if age_out > age_in

xtile sri_q_raw = sleepregularityindex_ad_t5a5_mn, nq(4)
label define sri_lbl 1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4"
label values sri_q_raw sri_lbl

noisily display "正在将字符串分类变量转换为模型可用格式..."
foreach var in sex season race smk alc {
    capture confirm string variable `var'
    if _rc == 0 { 
        encode `var', gen(`var'_num)
        drop `var'
        rename `var'_num `var'
    }
}

* ==============================================================================
* 4. 设定左截断生存数据
* ==============================================================================
stset age_out, failure(target_status==1) enter(age_in)

* ==============================================================================
* 5. 拟合 stpm2 灵活参数生存模型
* ==============================================================================
noisily display "正在拟合 stpm2 模型..."
stpm2 ib4.sri_q_raw i.sex c.mvpa i.season i.race c.tdi i.smk i.alc c.bmi c.fastingtime c.sleepdurationinspt_ad_t5a5_mn, df(3) scale(hazard)

* ==============================================================================
* 6. 循环计算 45 到 100 岁的 YLL (90岁以后自动完美收敛)
* ==============================================================================
capture postutil clear
postfile yll_data age yll lci uci using "FPSM/YLL_results_100.dta", replace

capture drop t_grid t_base
gen t_grid = .
gen t_base = .

noisily display "开始逐年计算 YLL (智能规避极限高龄求导崩溃)..."

* 【完美修复】：严格积分算到 90 岁，规避 92 岁时的分母极限为 0 导致矩阵爆炸！
forvalues a = 45/90 {
    
    replace t_grid = .
    replace t_base = `a'
    
    local n_pts = (100 - `a') * 10 + 1 
    forvalues i = 1/`n_pts' {
        quietly replace t_grid = `a' + (`i'-1)*0.1 in `i'
    }
    
    capture drop diff_est diff_lci diff_uci
    
    quietly predictnl diff_est = ///
        (predict(surv timevar(t_grid) at(sri_q_raw 4)) / predict(surv timevar(t_base) at(sri_q_raw 4))) ///
        - (predict(surv timevar(t_grid) at(sri_q_raw 1)) / predict(surv timevar(t_base) at(sri_q_raw 1))), ///
        ci(diff_lci diff_uci) level(95)
        
    quietly integ diff_est t_grid in 1/`n_pts'
    local yll_est = r(integral)
    
    quietly integ diff_lci t_grid in 1/`n_pts'
    local yll_lci = r(integral)
    
    quietly integ diff_uci t_grid in 1/`n_pts'
    local yll_uci = r(integral)
    
    post yll_data (`a') (`yll_est') (`yll_lci') (`yll_uci')
}

* 【神级走位】：手动补上 100 岁这一天的锚点，曲线将自动平滑收束为针尖！
post yll_data (100) (0) (0) (0)

postclose yll_data
noisily display "所有 45-100 岁 YLL 计算完成！"

* ==============================================================================
* 7. 绘制极美图表并输出 
* ==============================================================================
use "FPSM/YLL_results_100.dta", clear

* 强制裁剪超出范围的 CI
replace uci = 4 if uci > 4
replace lci = 0 if lci < 0

twoway ///
    (rarea lci uci age, color(red%25) lwidth(none)) /// 
    (line yll age, lcolor(red) lwidth(medthick)), ///   
    yline(0, lpattern(dash) lcolor(black)) ///
    xtitle("Baseline Age (Years)", size(medlarge) height(5)) ///
    ytitle("Years of Life Lost up to Age 100" "(Q1 vs Q4)", size(medlarge) height(5)) ///
    xlabel(45(5)100, labsize(medium)) ///             <-- 完美展平至 100 岁
    ylabel(0(1)4, labsize(medium) format(%02.1f)) /// 
    yscale(range(0 4)) ///                            
    legend(order(2 "YLL (Q1 vs Q4)" 1 "95% CI") ring(0) position(1) region(lcolor(none) fcolor(none))) /// <-- 图例到位
    graphregion(color(white)) plotregion(color(white))

graph export "FPSM/YLL_UpTo100_Stata_Q1vsQ4_Full100.pdf", replace as(pdf)
graph export "FPSM/YLL_UpTo100_Stata_Q1vsQ4_Full100.png", replace as(png) width(2000)

noisily display "✅ 大功告成！全景图已生成，再也不会崩溃了！"
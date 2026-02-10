import numpy as np
import pandas as pd
import os

# 결과 파일 경로 설정
input_file = r'C:\Users\daybr\Git_Work\CMOST_experiment\python\Results\CMOST13_02072016_Results.npz'
output_file = input_file.replace('.npz', '.csv')

if os.path.exists(input_file):
    print(f"Loading {input_file}...")
    data = np.load(input_file, allow_pickle=True)
    
    # 주요 데이터 추출
    if 'Variable' in data and 'Var_Legend' in data:
        variables = data['Variable']
        legends = data['Var_Legend']
        
        # 데이터 길이 맞추기 (오류 방지)
        min_len = min(len(variables), len(legends))
        
        # DataFrame 생성
        df = pd.DataFrame({
            'Description': legends[:min_len],
            'Value': variables[:min_len]
        })
        
        # CSV로 저장 (한글 깨짐 방지 utf-8-sig)
        df.to_csv(output_file, index=False, encoding='utf-8-sig')
        print(f"Succeessfully saved to {output_file}")
        
        # 내용 미리보기
        print("\n--- Preview ---")
        print(df.head(10))
    else:
        print("Error: 'Variable' or 'Var_Legend' not found in the file.")
else:
    print(f"File not found: {input_file}")
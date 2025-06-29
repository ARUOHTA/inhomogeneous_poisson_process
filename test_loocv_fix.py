#!/usr/bin/env python3
"""
修正されたLOOCVシステムの小規模テスト
"""

import sys
from pathlib import Path

# プロジェクトルートをPythonパスに追加
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

import numpy as np
from bayesian_statistics.model3_config import Model3Config, Model3Pipeline
from bayesian_statistics.model3_loocv import LOOCVConfig, LOOCVEvaluator


def main():
    print("=== 修正されたLOOCVシステムのテスト ===")
    
    # 設定
    data_dir = "/home/ohta/dev/bayesian_statistics/data/"
    
    model_config = Model3Config(data_dir=data_dir)
    loocv_config = LOOCVConfig(
        n_trials=2,  # 2試行でテスト
        verbose=True,
        random_seed=42
    )
    
    try:
        # 前処理
        print("前処理を実行中...")
        pipeline = Model3Pipeline(model_config)
        preprocessor = pipeline.run_preprocessing()
        
        # LOOCV評価器の作成
        print("LOOCV評価器を初期化中...")
        evaluator = LOOCVEvaluator(
            preprocessor=preprocessor,
            model_config=model_config,
            loocv_config=loocv_config
        )
        
        # 利用可能な遺跡を確認
        test_period = 0
        available_sites = []
        for site_id in evaluator.site_ids[:20]:  # 最初の20個から探す
            if site_id in evaluator.observed_compositions[test_period]:
                comp = evaluator.observed_compositions[test_period][site_id]
                if not np.allclose(comp, 1.0 / len(comp)):
                    available_sites.append(site_id)
        
        if len(available_sites) == 0:
            print("テスト可能な遺跡が見つかりませんでした")
            return
        
        test_site_id = available_sites[0]
        print(f"テスト遺跡: {test_site_id}")
        print(f"観測構成比: {evaluator.observed_compositions[test_period][test_site_id]}")
        
        # 複数試行を実行
        print(f"\n{loocv_config.n_trials}試行のLOOCVを実行中...")
        period_result = evaluator.run_period_evaluation(test_period, n_trials=loocv_config.n_trials)
        
        print(f"\n=== 結果サマリー ===")
        stats = period_result.get('summary_statistics', {})
        print(f"実行試行数: {period_result.get('n_trials', 0)}")
        print(f"成功試行数: {stats.get('n_successful_trials', 0)}")
        print(f"成功率: {stats.get('success_rate', 0.0):.3f}")
        
        if stats.get('n_successful_trials', 0) > 0:
            print(f"平均Aitchison距離: {stats.get('mean_aitchison_distance', np.nan):.4f} ± {stats.get('std_aitchison_distance', np.nan):.4f}")
            print(f"平均Total Variation: {stats.get('mean_total_variation', np.nan):.4f} ± {stats.get('std_total_variation', np.nan):.4f}")
            
            # 最初の3試行の詳細を表示
            print(f"\n=== 詳細結果（最初の3試行） ===")
            for i, trial in enumerate(period_result['trial_results'][:3]):
                if trial['success']:
                    print(f"\n試行 {i+1} (遺跡ID: {trial['test_site_id']}):")
                    print(f"  Aitchison距離: {trial['aitchison_distance']:.4f}")
                    print(f"  Total Variation: {trial['total_variation']:.4f}")
                    print(f"  観測: {trial['observed_composition']}")
                    print(f"  予測: {trial['predicted_composition']}")
                else:
                    print(f"\n試行 {i+1}: 失敗 - {trial['error']}")
        
    except Exception as e:
        print(f"エラーが発生しました: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()

## ファイル

scripts/ : ベンチマーク用スクリプト
docs/ : outline.txtが入っている

blast.cc : blast_SemiGappedAlign
simdblast.cc : blast_SIMDSemiGappedAlign
rognes.cc : static banded with horizontal packed band
diag.cc : static banded with anti-diagonal packed band
ddiag.cc : dynamic banded with anti-diagonal packed band

## コンパイル

### 計算時間のベンチ

	$ g++ -O3 -Wall main.cc diag.cc ddiag.cc blast.cc simdblast.cc rognes.cc util.cc -msse4.1 -DBENCH -DBW=32 -DALL
	$ ./a.out

* -O1以上で最適化しないと _mm_extract_epi16 が即値しか受け取れないという
エラーが出るので、必ず-O1以上をつける
* -DBENCH をつけると時間を計測するマクロが埋め込まれる
* -DBW=xx でバンド幅を指定
* -DALL で全ての (diag, rognes, simdblast) ベンチマークを実行


### パラメータを振るベンチ

スコア行列を振るベンチはサフィックスなし、ギャップ長さを振るベンチは
_gapのサフィックスがついたスクリプトをつかう。

#### コンパイル

	$ for b in 8 16 32 48 56 64; do for f in blast diag ddiag simdblast rognes; do clang++ -O3 -Wall $f.cc util.cc -msse4.1 -DMAIN -DBW=$b -o ./bin/$f-$b; done; done

#### 1スレッドで走らせるとき

	$ python scripts/generate_params.py > params.tsv
	$ python scripts/evaluate.py params.tsv out.tsv
	$ python scripts/aggregate.py out.tsv result.txt

でシミュレーションが走る。

	$ python
	>>> from scripts.util import *
	>>> (l, a) = load_result('result.txt')

で結果が l, a に numpy.array で読み込まれる。l[:, :, 3, 3, 3]などで
テーブルを参照。それぞれのフィールドはscripts/params.pyのparams_listを
参照。

#### SGEで走らせるとき

	$ scripts/sim.sh

でジョブを投入する。

	$ scripts/aggr.sh

で結果をまとめる。1スレッドで走らせたときと同じ result.txt が生成される。


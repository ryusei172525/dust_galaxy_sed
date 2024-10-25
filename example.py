import asyncio
from galaxy_sed.mcmc import MCMCOptimizer
from galaxy_sed.diagram import delete_past_data

galaxy_age = 5  # 4: 500Myr, 5: 600Myr

async def main():
    delete_past_data()
    result = await MCMCOptimizer.run_mcmc(galaxy_age)  # 非同期関数内でawaitを使用
    print("計算結果:", result)

# イベントループを実行
if __name__ == "__main__":
    asyncio.run(main())  # main関数を非同期に実行
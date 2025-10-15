import calendar
import datetime

import jpholiday


def count_holidays_in_month(year: int, month: int) -> int:
    """
    指定された年と月の中で、土日と日本の祝日を合わせた休日の数を返す。
    ただし、祝日が土日と重なった場合は重複しないように扱う。
    """
    holiday_dates = set()

    # 月の最終日を取得
    _, last_day = calendar.monthrange(year, month)

    for day in range(1, last_day + 1):
        date = datetime.date(year, month, day)

        # 土日
        if date.weekday() >= 5:  # 5: 土曜, 6: 日曜
            holiday_dates.add(date)

        # 祝日
        if jpholiday.is_holiday(date):
            holiday_dates.add(date)

    return len(holiday_dates)

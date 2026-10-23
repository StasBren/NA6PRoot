#!/usr/bin/env python3
"""Read-only localisation of verified 100k-parent six-geometry MWPC HTML.

Localises rendered pages (not physics) without rerunning Geant4, Heed or
Garfield. Preserves the embedded machine-readable numeric JSON byte-for-byte.
Creates *_en.html siblings, leaving Russian HTML and frozen merged JSON alone.
Every generated English page is audited for untranslated Cyrillic and all
six fixed analysis configurations retain the same payload and navigation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from pathlib import Path

from run_pi_geant4_campaign import REPO

DEFAULT=REPO/"test_runs/mwpc_stage_c/pi_six_geometry_2m_timing_v1/full/pages"
SOURCES=("jpsi40-bore","phi40-bore","omega40-bore")
CYRILLIC=re.compile(r"[А-Яа-яЁё]")
DATA_TAG=re.compile(r'(<script type="application/json" id="(?:payload|data)">)(.*?)(</script>)',re.S)
TRANSLATIONS=[["NA60+ — PI · XY-карта, время дрейфа и углы · 6 геометрий","NA60+ · MWPC | XY map, electron drift time and arrival angles | Six geometries"],["NA60+ · PI · NATIVE GEANT4 + HEED + GARFIELD · ШЕСТЬ ГЕОМЕТРИЙ","NA60+ · MWPC · NATIVE GEANT4 + HEED + GARFIELD · SIX GEOMETRIES"],["Время дрейфа, аноды и углы — сравнение 3 / 3.5 / 4 мм","Electron drift, anodes and arrival angles — 3 / 3.5 / 4 mm comparison"],["Физическая геометрия камер сохраняется; шаг анодов и заземлённые промежуточные проволоки диаметром 75 мкм переключаются без новой симуляции. Трёхтысячный первый чанк J/ψ / ϕ / ω на источник.","The physical chamber layout is unchanged. Compare anode pitches with or without grounded 75 µm field wires using the same underlying muon sample."],["Трёхтысячный первый чанк J/ψ / ϕ / ω на источник.","First 3,000 generated J/ψ / ϕ / ω parent events per source."],["сгенерированных родительских событий на источник; объединены проверенные независимые Geant4-чанки.","generated parent events per source; independently validated Geant4 chunks combined."],["истинные native Geant4 газовые пересечения из независимых чанков","true native Geant4 gas crossings from independently validated chunks"],["сгенерированных событий на канал","generated events per source"],["Оригинальные отчёты на 3 000 событий сохранены отдельно.","The previous 3,000-parent-event reports are preserved separately."],["← Навигация по источникам и станциям","← Source / station index"],["Исходный V2, только 4 мм","Earlier V2 (Russian; 4 mm only)"],["Детальные спектры отдельных анодов","Detailed per-anode spectra"],["Версия:","Version:"],["Выбор камеры, угла и проволок","Chamber, angle and wire controls"],["Анодный шаг и полевая проволока","Anode pitch and field wire"],["Угол входа мюона","Muon incidence angle"],["θu — поперёк анодов","θu — across anode wires"],["θw — вдоль анодов","θw — along anode wires"],["Интервал угла (только с данными)","Angle bin (only bins with data)"],["Конкретный анод выбранной камеры","Physical anode within selected chamber"],["Вертикальная шкала времени","Drift-time vertical scale"],["Линейная","Linear"],["Логарифмическая","Logarithmic"],["Время: сравнить конфигурации","Drift-time configuration comparison"],["Один вариант","One configuration"],["Все 6 вариантов","All six configurations"],["Угол прихода электрона","Electron arrival angle"],["К нормали проволоки","Relative to wire-surface normal"],["Положение попадания на проволоке φ","Wire-surface impact position angle φ"],["Выбор камеры — кликом на XY-карту.","Select a chamber by clicking its rectangle on the XY map."],["Угол входа мюона сохраняет знак. Временной PDF относится к первичным электронам, собранным анодами с известным временем, включая хвост ≥300 нс. Ниже — не аналоговый сигнал и не время срабатывания ASIC.","Muon entrance angles retain their signs. The drift-time PDF is conditional on primary electrons collected at anodes with known drift times, including the ≥300 ns overflow bin. This is neither an analogue induced-current waveform nor the ASIC trigger time."],["В этом наборе отсутствуют","Not available in this dataset:"],["распределения рождения электронов по времени и P(N≥2): старый V2 считал их для другой, единственной конфигурации. Они не подменяются модельными величинами.","electron birth-time distributions and event-level P(N≥2) from the earlier V2 geometry. They are not replaced by synthetic estimates."],["1. Выбор камеры на XY-плане станции","1. Physical chamber selection — station XY map"],["Кликните на прямоугольник DET. Цвет — доля первичных электронов, достигших анода за первые 100 нс, относительно всех первичных электронов из газовых сегментов этой камеры и выбранного углового бина.","Click a DET rectangle. The colour encodes primary electrons collected at an anode within 100 ns, divided by all primary electrons generated in the selected chamber and muon-angle bin."],["Это не P(N≥2).","This is NOT P(N≥2)."],["Серый — недостаточно данных. Оси XY в сантиметрах.","Grey indicates no primary-electron denominator. X and Y are in centimetres."],["Все камеры станции","All station chambers"],["Ожидание загрузки интерактивной XY-карты…","Loading interactive XY chamber map…"],["Карта физических камер с выбором DET и осями X,Y, сантиметры","Selectable physical chamber layout, X and Y axes in centimetres"],["2. Время дрейфа — все камеры","2. Electron drift time — all chambers"],["Время от рождения электрона до сбора на аноде. Один прямоугольник — 10 нс.","Elapsed time from primary-electron birth to collection at an anode. Each bin spans 10 ns."],["Время дрейфа от рождения электрона, проценты по 10 нс","Drift-time distribution of primary electrons, percent per 10 ns"],["Нормировка каждого варианта: все собранные первичные электроны с известным временем, включая хвост ≥300 нс. Отбор по углу входа мюона действует и для времени, и для карты камер.","Each configuration is normalised to all anode-collected primary electrons with known drift time, including the ≥300 ns overflow. Muon incidence-angle selection applies to both the time histogram and the chamber map."],["3. Как электроны подходят к аноду?","3. Electron arrival geometry at the anode"],["Угол между последним микроскопическим шагом Garfield и внутренней нормалью к аноду. Не угол входа мюона.","Angle between Garfield's final microscopic electron step and the inward anode-surface normal, not the incoming muon angle."],["Распределение углов прихода электронов","Electron arrival-angle distribution"],["Бины по 10°. Выбор станции, камеры, анода, поля и знакового угла мюона синхронизован со временем. Для отдельного анода доступна только интеграция по всем углам входа мюона.","10° bins. Station, chamber, anode, field configuration and signed muon angle follow the time-histogram selection. Per-anode distributions are available only integrated over muon entrance angles."],["4. Шесть вариантов: сбор, Mean и centered RMS","4. Six geometries: anode collection, mean and centred RMS"],["Те же исходные Geant4/Heed-события и 6 разных Garfield LUT. Средние по","The same underlying Geant4/Heed events are folded through six different Garfield LUTs. Means are calculated for"],["собранным","collected"],["электронам — не весь исходный поток; эффективность может меняться.","electrons, not all primary electrons; collection efficiency may change."],["Внимание:","Important:"],["один и тот же номер физического анода при разных шагах относится к разным координатам внутри камеры; для сопоставления одинаковой области лучше выбирать всю камеру.","the same physical anode index corresponds to different chamber coordinates at different pitches; select the whole chamber for matched-area comparisons."],["Среднее время дрейфа шести геометрий","Mean drift time across six geometries"],["p / проволока","Pitch / field wire"],["Электронов","Electrons"],["Доля сбора*","Collection fraction*"],["Mean, нс","Mean, ns"],["RMS, нс","RMS, ns"],["*Доля сбора по всем первичным электронам определяется на уровне камеры/станции. Для выбранного физического анода знаменатель первичных электронов, принадлежащих именно аноду, не определён. Числа — без bootstrap-ошибок; отдельный event-level отчёт будет подготовлен.","*The fraction of all primary electrons collected at anodes is defined at chamber/station level. A primary-electron denominator uniquely attributable to one physical anode is not defined. Parent-event bootstrap uncertainties have not yet been calculated."],["Данные: истинные native Geant4 газовые пересечения первого чанка (3 000 сгенерированных событий на канал), Heed-ионизация, B=0 Garfield LUT 3+3 мм, 30 мкм аноды; шаг 3 / 3.5 / 4 мм, без и с заземлёнными 75 мкм промежуточными проволоками.","Data: true native Geant4 gas crossings (3,000 generated parents per source), Heed ionisation and B=0 Garfield LUTs; symmetric 3+3 mm gap, 30 µm anodes, 3 / 3.5 / 4 mm pitches, with or without grounded 75 µm field wires."],["Данные:","Data:"],["Heed-ионизация, B=0 Garfield LUT 3+3 мм, 30 мкм аноды; шаг 3 / 3.5 / 4 мм, без и с заземлёнными 75 мкм промежуточными проволоками.","Heed ionisation, B=0 Garfield LUTs, symmetric 3+3 mm gas gap and 30 µm anodes; pitches 3 / 3.5 / 4 mm, with or without grounded 75 µm field wires."],["18 детальных страниц и исходные 496 809 групп не изменяются.","The 18 detailed pages and the original 496,809 groups are unchanged."],["Не являются расчетом ASIC/лавинного сигнала и пока не имеют parent-event статистических ошибок по детектору и углу.","This is not an avalanche-current or ASIC-response simulation. Parent-event statistical errors by chamber and angle are not yet provided."],["Шесть геометрий","Six geometries"],["Источник","Source"],["Станция","Station"],["Все углы (","All angles ("],[" хитов)"," gas hits)"],["Все аноды","All anodes"],["Анод ","Anode "],["анод ","anode "],["Все камеры","All chambers"],["без field wire","no field wire"],["field wire 75 мкм","75 µm field wire"],["нет данных","no data"],["P(e⁻≤100нс)=","P(e⁻ collected in 100 ns)="],["X, см →","X, cm →"],["Y, см ↑","Y, cm ↑"],[" · газовых хитов "," · gas hits "],[" · анод "," · anode "],["Цветовая шкала: 0–","Colour scale: 0–"],[" P(первичный электрон собран за первые 100 нс). Серый: нет первичных e⁻. "," P(primary electron collected within 100 ns). Grey: no primary electrons. "],["Жёлтая рамка: выбранная камера. НЕ P(N≥2).","Yellow outline: selected chamber. NOT P(N≥2)."],["Газовых хитов","Muon gas crossings"],["Первичных электронов","Primary electrons"],["Собрано","Collected"],["Время известно","Known drift time"],["Время неизвестно","Unknown drift time"],["Потерь на краях","Finite-edge losses"],["Mean, нс","Mean, ns"],["Centered RMS, нс","Centred RMS, ns"],["≥300 нс","≥300 ns"],["Время дрейфа от рождения до анода, нс","Electron birth-to-anode drift time, ns"],["Доля e⁻ / 10 нс, % (логарифмическая шкала)","Electrons per 10 ns, % (logarithmic axis)"],["Доля e⁻ в интервале 10 нс, %","Electrons per 10 ns, %"],["нс: ","ns: "],[" электронов, "," electrons, "],["100 нс","100 ns"],["Нет собранных электронов с известным временем","No collected electrons with known drift time"],["2. Время дрейфа — ","2. Electron drift time — "],["все камеры","all chambers"],[", анод ",", anode "],["Угол подхода электрона к внутренней нормали поверхности анода; направление последнего микроскопического шага Garfield.","Electron approach to the inward anode-surface normal; direction of the last microscopic Garfield step."],["Геометрический угол положения попадания на поверхность анода φ (со знаком), а НЕ угол направления электрона.","Signed geometric impact-position angle φ on the anode surface, NOT the direction of electron motion."],["Угол к нормали анода, градусы","Angle relative to anode normal, degrees"],["Угол положения попадания φ, градусы","Impact-position angle φ, degrees"],["Доля электронов в интервале 10°, %","Electron fraction per 10°, %"],["Для выбранного фильтра нет известных углов","No known arrival angles for the selected filter"],["Конфигурации p, мм / полевая проволока","Anode pitch p, mm / field wire"],["Mean t дрейфа, нс","Mean drift time, ns"],["Нет подходящих данных для выбора; измените фильтр.","No matching data; change the selection."],["Ошибка интерфейса: ","Interface error: "],[". Резервный статический спектр ниже остаётся доступным.",". The static fallback spectrum remains available below."],["NA60+ — время дрейфа для разных анодов","NA60+ — per-anode electron drift-time distributions"],["NA60+ · Дрейф электронов к анодным проволокам","NA60+ · Primary electron drift to anode wires"],["Версия интерфейса:","Interface version:"],["Реальные локальные газовые пересечения Geant4 → общие для шести вариантов","True local Geant4 gas crossings → shared primary-electron births for all six"],[" первичные электроны Heed → Garfield LUT. Время отсчитывается от рождения"," configurations (Heed) → Garfield LUT. Time is measured from electron birth"],["электрона, а не от пролёта мюона или порога ASIC.","rather than muon arrival or ASIC threshold crossing."],["Камера","Chamber"],["Физический анод","Physical anode"],["Шаг и полевая проволока","Anode pitch and field wire"],["Интервал угла","Angle bin"],["Наблюдаемая величина","Observable"],["Время дрейфа [нс]","Drift time [ns]"],["Угол прихода относительно нормали [°]","Electron arrival angle relative to normal [°]"],["Угол положения попадания на аноде [°]","Anode-surface impact-position angle [°]"],["Вертикальная шкала","Vertical scale"],["Отображение","Plot mode"],["Гистограмма выбранной конфигурации","Histogram of selected configuration"],["Наложение всех шести конфигураций","Overlay all six configurations"],["Ожидание JavaScript. Если график не появился, статический спектр ниже доступен без скриптов.","JavaScript is loading. If the chart does not render, a static reference spectrum is shown below."],["Доля электронов в интервале 10 нс, %","Electron fraction per 10 ns, %"],["Спектр времени дрейфа с шагом 10 наносекунд","Drift-time spectrum with 10 ns bins"],["Временная нормировка включает ≥300 нс; они вынесены","The drift-time normalisation includes ≥300 ns electrons, recorded"],["в отдельный счётчик. Угловые гистограммы используют электроны с известным","in a separate overflow count. Angular spectra use electrons with known"],["углом, что не обязательно совпадает с выборкой известных времен.","angles, which need not be the same population as those with known drift times."],["Штриховая линия 100 нс показывается только на временном графике.","The 100 ns reference line appears on drift-time charts only."],["Угол попадания φ — геометрический угол на поверхности проволоки, а","The impact-position angle φ describes position on the wire surface, while"],["угол к нормали — направление последнего микроскопического шага Garfield.","the incidence angle describes Garfield's final electron step direction."],["Пустые группы не заменяются данными других камер.","Missing groups are not filled from other chambers."],["Границы физической интерпретации","Limits of physical interpretation"],["В режиме «анод» число газовых сегментов означает сегменты с ≥1","In anode mode, the gas-segment count means segments contributing at least one"],["собранным электроном на выбранном аноде; это не количество входящих","collected electron to that anode; it is not the number of muons entering"],["мюонов в его геометрическую область. Статистические ошибки по родительским","its geometric area. Parent-event statistical uncertainties have"],["Geant4-событиям ещё не вычислены. Mean и centered RMS ниже — описательные","not yet been calculated. Mean and centred RMS below are descriptive"],["оценки данного первого чанка, не приборная временная разрешающая способность.","statistics, not instrumental time resolution."],["оценки выборки из ","estimates for a sample of "],[" сгенерированных родительских событий на источник"," generated parent events per source"],["Все источники","All sources"],["Без отбора","No angle selection"],["θu (со знаком)","θu (signed)"],["θw (со знаком)","θw (signed)"],["Все углы","All angles"],["DET ","DET "],["анод ","anode "],["Газовых сегментов","Muon gas crossings"],["Собрано на анод","Collected at anode"],["Неизвестное время","Unknown drift time"],["Потери на конечных краях","Finite-edge losses"],[" — выбранная группа отсутствует"," — selected group has no data"],["Время дрейфа от рождения электрона, нс","Electron birth-to-anode drift time, ns"],["Угол подхода электрона к нормали, градусы","Electron approach angle relative to normal, degrees"],["Геометрический угол попадания φ, градусы","Geometric impact-position angle φ, degrees"],["Доля электронов в интервале 10°, %","Electron fraction per 10°, %"],["Физическая ширина газового объёма сохраняется; анодов на станциях: ","Physical gas width is unchanged; anodes per station: "],["физическая ширина газового объёма сохраняется; анодов на станциях: ","physical gas width unchanged; anodes per station: "],["Недостаточно элементов отрисовано: stats ","Insufficient rendered elements: stats "],["Ошибка построения графика: ","Chart rendering error: "],[". Исходные физические данные сохранены; сообщи текст ошибки.",". The underlying physics data are unchanged; report the error text."],["Резервная гистограмма — построена Python без JavaScript","Static fallback histogram — generated without JavaScript"],["Если интерактивный график остаётся пустым, этот результат ","If the interactive chart remains blank, this static plot "],["отображает действительные данные для выбранной станции и ","shows actual data for the selected station and "],["конфигурации p=3 мм без полевой проволоки.","configuration p=3 mm without a field wire."],["Конфигурация","Configuration"],["Известное время","Known drift time"],["Гистограмма времени дрейфа, резервный график","Static drift-time histogram"],["Время дрейфа электрона, нс","Electron drift time, ns"],["Упрощённый обзор","Simplified overview"],["Первоначальный подробный дашборд","Original detailed dashboard (Russian)"],["NA60+ · PI · шесть конфигураций проволок","NA60+ · MWPC · Six wire geometries"],["NA60+ · время дрейфа · 6 геометрий проволок","NA60+ · Electron drift time · Six wire geometries"],["Интерфейс как в исходном pi_centered_timing_dashboard_v2:","The interactive interface follows the previous pi_centered_timing_dashboard_v2:"],["кликабельная XY-карта физических камер слева и время дрейфа справа.","physical XY chamber map on the left and electron drift time on the right."],["Шаг p=3, 3.5, 4 мм; без и с заземлённой полевой проволокой","Anode pitch p=3, 3.5, 4 mm, with or without a grounded field wire"],["диаметром 75 мкм. Ниже — углы прихода электронов и сравнительная таблица.","of 75 µm diameter. Electron arrival angles and configuration comparisons are below."],["Источник и станцию можно менять на каждой странице.","Source and station can be changed on each page."],["интерактивная XY-карта камер","interactive XY chamber map"],["Данные: реальный нативный первый чанк Geant4/Heed + шесть Garfield","Data: first verified native Geant4/Heed chunk + six Garfield"],["сгенерированных родительских событий на источник; Geant4/Heed + шесть Garfield","generated parents per source; Geant4/Heed + six Garfield"],["LUT; не пороговая электроника. Базовые результаты не изменяются.","LUTs; no discriminator simulation. Source results are unchanged."],["ШЕСТЬ ГЕОМЕТРИЙ","SIX GEOMETRIES"],["первичные электроны Heed → Garfield LUT. Время отсчитывается от рождения","primary electrons from Heed → Garfield LUT. Time is measured from the birth of the"],["Mean [нс]","Mean [ns]"],["Centered RMS [нс]","Centred RMS [ns]"],["≥ 300 нс","≥ 300 ns"],["Газовых сегментов","Muon gas crossings"],["Электронов","Electrons"],["Отображение","Display mode"],["Физическая геометрия камер сохраняется; шаг анодов и заземлённые промежуточные проволоки диаметром 75 мкм переключаются без новой симуляции.","The physical chamber layout is unchanged; the anode pitch and grounded 75 µm intermediate field wires are compared using the same muon sample."],["(по ","("],[", не приборная временная разрешающая способность.",", not the instrumental timing resolution."]]

def sha(data):
    return hashlib.sha256(data.encode("utf-8")).hexdigest()

def translate_page(src,dst,kind,*,write=True):
    if not src.is_file():
        raise FileNotFoundError(src)
    html=src.read_text(encoding="utf-8")
    embedded=[]
    def hide(match):
        i=len(embedded)
        embedded.append(match.group(0))
        return f"__ENGLISH_NUMERIC_PAYLOAD_{i}__"
    stripped=DATA_TAG.sub(hide,html)
    if kind in ("layout","detail") and len(embedded)!=1:
        raise ValueError(f"Expected EXACTLY one frozen analysis JSON payload: {src}")
    if kind in ("layout","detail"):
        # The per-anode pages can embed tens of megabytes of group JSON.
        # Check frozen exposure directly in the metadata text without
        # deserialising every electron bin into a second Python object tree.
        data_text=DATA_TAG.search(html).group(2)
        exposure=re.search(
            r'"n_generated_parent_events_per_source"\s*:\s*(\d+)',
            data_text[:131072])
        nparents=int(exposure.group(1)) if exposure else None
        if nparents!=100_000:
            raise ValueError(
                f"Refusing English 100k label: true parent exposure is {nparents}"
                f" in {src}")
    if kind=="portal" and embedded:
        raise ValueError("Portal unexpectedly embeds JSON analysis")
    if kind=="layout":
        stripped=stripped.replace("../pi_pitch_fieldwire_timing_dashboard_v3.html",
                                  "../pi_pitch_fieldwire_timing_dashboard_v3_en.html")
        stripped=stripped.replace('MS"+S.st+".html"','MS"+S.st+"_en.html"')
        stripped=stripped.replace('MS"+S.st+"_layout_v3.html"',
                                  'MS"+S.st+"_layout_v3_en.html"')
        # Historical V2 was Russian; its link remains, clearly labelled.
    elif kind=="portal":
        stripped=re.sub(r'(<a class="card" href="[^"]+/MS\d+_layout_v3)\.html"',
                        r'\1_en.html"',stripped)
        # The Russian original points to two files absent in 100k/pages.
        stripped=re.sub(r'<p><a href="index.html">.*?</a></p>',"",stripped,flags=re.S)
    # Detailed-page JS is unchanged except its written English labels and nav.
    # Long phrases first, so replacing a short common noun cannot split them.
    for before,after in sorted(TRANSLATIONS,key=lambda x:len(x[0]),reverse=True):
        stripped=stripped.replace(before,after)
    stripped=stripped.replace('lang="ru"','lang="en"')
    stripped=stripped.replace('"ru-RU"','"en-US"')
    # The portal's V2 links may legitimately lead to earlier Russian pages,
    # but the localized pages themselves must be fully English.
    issues=[]
    for m in list(CYRILLIC.finditer(stripped))[:6]:
        issues.append(repr(stripped[max(0,m.start()-55):m.start()+85]))
    if issues:
        raise ValueError(f"Untranslated Cyrillic in {src}:\n"+"\n".join(issues))
    for i,payload in enumerate(embedded):
        stripped=stripped.replace(f"__ENGLISH_NUMERIC_PAYLOAD_{i}__",payload)
    if kind in ("layout","detail") and DATA_TAG.findall(stripped)!=DATA_TAG.findall(html):
        raise AssertionError("Physics JSON changed during translation")
    if kind=="layout":
        for marker in ("id=\"map\"","id=\"time\"","id=\"arrival\"",
                       "id=\"comparison\"","timeScale"):
            if marker not in stripped:
                raise ValueError(f"Missing original V3 interactive element {marker}")
    if kind=="detail" and 'id="yscale"' not in stripped:
        raise ValueError("Lost detailed time-spectrum log-Y selector")
    if write:
        dst.write_text(stripped,encoding="utf-8")
    return {"path":str(dst),"source_sha256":sha(html),
            "english_sha256":sha(stripped),"bytes":len(stripped.encode("utf-8")),
            "preserved_json_payloads":len(embedded)}

def build(base,*,dry_run=False):
    base=base.resolve()
    if not base.is_dir():raise NotADirectoryError(base)
    produced=[]
    for source in SOURCES:
        for st in range(6):
            folder=base/source
            for suffix,kind in (("_layout_v3","layout"),("","detail")):
                stem=f"MS{st}{suffix}"
                produced.append(translate_page(
                    folder/f"{stem}.html",folder/f"{stem}_en.html",kind,
                    write=not dry_run))
    portal="pi_pitch_fieldwire_timing_dashboard_v3"
    produced.append(translate_page(
        base/f"{portal}.html",base/f"{portal}_en.html","portal",
        write=not dry_run))
    manifest={"version":"PI_NATIVE_100K_ENGLISH_HTML_V1",
              "translation_only":True,
              "physics_simulation_performed":False,
              "generated_parent_events_per_source":100_000,
              "english_pages":len(produced),"artifacts":produced}
    if not dry_run:
        file=base/"english_dashboard_manifest.json"
        file.write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",
                        encoding="utf-8")
    print(f"PASS: {len(produced)} English V3, detailed spectra and portal HTML "
          f"pages; numeric embedded JSON unchanged; "
          f"{'dry-run' if dry_run else 'files written'}")
    print("English portal:",base/f"{portal}_en.html")
    return manifest

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pages",type=Path,default=DEFAULT)
    p.add_argument("--dry-run",action="store_true")
    args=p.parse_args()
    build(args.pages,dry_run=args.dry_run)

if __name__=="__main__":
    main()

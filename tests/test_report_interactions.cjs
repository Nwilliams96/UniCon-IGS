/* Optional real-browser regression tests for the self-contained report.
 * Requires Node and Playwright. First call render_interaction_fixture() from
 * test_generate_pipeline_report.py, then set REPORT_FIXTURE_HTML to its output.
 * REPORT_BROWSER_EXECUTABLE may select an installed, compatible Chrome binary.
 * No sequencing tools, real study data, network services or user profiles are used.
 */
const {test, before, after} = require("node:test");
const assert = require("node:assert/strict");
const {readFileSync, mkdirSync} = require("node:fs");
const {pathToFileURL} = require("node:url");
const {join} = require("node:path");
const {chromium} = require("playwright");

const fixture = process.env.REPORT_FIXTURE_HTML;
assert.ok(fixture, "Set REPORT_FIXTURE_HTML to render_interaction_fixture()'s output");
let browser;

before(async () => {
  browser = await chromium.launch({
    headless: true,
    executablePath: process.env.REPORT_BROWSER_EXECUTABLE || undefined,
    args: ["--disable-crash-reporter"],
  });
});
after(async () => { if (browser) await browser.close(); });

async function withReport(run, options = {}) {
  const context = await browser.newContext({viewport: {width: 1280, height: 900}, ...options});
  const page = await context.newPage();
  const errors = [];
  const requests = [];
  page.on("pageerror", error => errors.push(String(error)));
  page.on("request", request => { if (/^https?:/.test(request.url())) requests.push(request.url()); });
  try {
    await page.goto(pathToFileURL(fixture).href);
    await run(page);
    assert.deepEqual(errors, [], "The report must have no browser exceptions");
    assert.deepEqual(requests, [], "The report must work offline");
  } finally {
    await context.close();
  }
}

async function sampleIds(page) {
  return page.locator(".taxonomy-sample-column").evaluateAll(columns => columns.map(column => column.dataset.sampleId));
}

async function selectFilterField(page, field) {
  await page.locator("#taxonomy-sample-filters").evaluate(panel => { panel.open = true; });
  await page.locator("#taxonomy-filter-field").selectOption(field);
}

function valueCheckbox(page, value) {
  return page.locator("#taxonomy-metadata-values input").filter({visible: true}).and(page.locator(`input[value=${JSON.stringify(value)}]`));
}

test("run summaries stay visible while both per-sample tables start collapsed", async () => withReport(async page => {
  for (const id of ["bbsplit-sample-details", "dada2-sample-details"]) {
    assert.equal(await page.locator(`#${id}`).evaluate(panel => panel.open), false);
    assert.equal(await page.locator(`#${id} table`).isVisible(), false);
    await page.locator(`#${id} summary`).click();
    assert.equal(await page.locator(`#${id} table`).isVisible(), true);
  }
  assert.equal(await page.locator("#bbsplit-sample-details tbody tr").count(), 5);
  assert.equal(await page.locator("#dada2-sample-details tbody tr").count(), 10);
  assert.match(await page.locator("#quality").innerText(), /10,000/);
}));

test("signed latitudes and decimal depths sort numerically with IDs preserved", async () => withReport(async page => {
  await page.locator("#taxonomy-plot-field").selectOption("Latitude");
  assert.deepEqual(await sampleIds(page), ["001", "002", "030", "010", "020"]);
  assert.equal(await page.locator(".taxonomy-x-axis-title").innerText(), "Latitude");
  assert.deepEqual(await page.locator(".taxonomy-sample-label").allTextContents(), ["-39.5", "-5.25", "12.25", "(blank)", "(blank)"]);
  await page.locator("#taxonomy-plot-field").selectOption("Depth");
  assert.deepEqual(await sampleIds(page), ["030", "002", "001", "010", "020"]);
  assert.equal(await page.locator(".taxonomy-x-axis-title").innerText(), "Depth");
  assert.deepEqual(await page.locator(".taxonomy-sample-label").allTextContents(), ["2.05", "2.5", "100.5", "(blank)", "(blank)"]);
  assert.match(await page.locator(".taxonomy-sample-label").first().getAttribute("title"), /^030\n/);
  assert.match(await page.locator(".taxonomy-sample-label").first().getAttribute("title"), /Temperature: 20.5/);
  await page.locator("#taxonomy-plot-field").selectOption("SampleID");
  assert.equal(await page.locator(".taxonomy-x-axis-title").innerText(), "SampleID");
  assert.deepEqual(await page.locator(".taxonomy-sample-label").allTextContents(), ["001", "002", "010", "020", "030"]);
}));

test("categorical and arbitrary metadata label the axis without changing sample identity", async () => withReport(async page => {
  await page.locator("#taxonomy-plot-field").selectOption("Condition");
  assert.equal(await page.locator(".taxonomy-x-axis-title").innerText(), "Condition");
  assert.deepEqual(await page.locator(".taxonomy-sample-label").allTextContents(), ["Environmental", "Environmental", "Environmental", "Mock", "(blank)"]);
  await page.locator("#taxonomy-plot-field").selectOption("Temperature");
  assert.equal(await page.locator(".taxonomy-x-axis-title").innerText(), "Temperature");
  assert.deepEqual(await page.locator(".taxonomy-sample-label").allTextContents(), ["20.5", "20.5", "21", "(blank)", "(blank)"]);
  assert.deepEqual(await sampleIds(page), ["001", "030", "002", "010", "020"]);
}));

test("translucent backgrounds do not fade plot contents and dropdown arrows are inset", async () => withReport(async page => {
  const style = await page.evaluate(() => {
    const control = document.querySelector(".explorer-control");
    const chart = document.querySelector(".frozen-chart");
    const bar = chart.querySelector("rect");
    return {
      arrowRight: getComputedStyle(control, "::after").right,
      arrowPointer: getComputedStyle(control, "::after").pointerEvents,
      background: getComputedStyle(chart).backgroundColor,
      chartOpacity: getComputedStyle(chart).opacity,
      barOpacity: getComputedStyle(bar).opacity,
      fill: bar.getAttribute("fill"), stroke: bar.getAttribute("stroke"),
    };
  });
  assert.equal(style.arrowRight, "14px");
  assert.equal(style.arrowPointer, "none");
  assert.equal(style.background, "rgba(255, 255, 255, 0.6)");
  assert.equal(style.chartOpacity, "1");
  assert.equal(style.barOpacity, "1");
  assert.equal(style.fill, "#24543d");
  assert.equal(style.stroke, "#111111");
}));

test("metadata exclusions are independent of sorting and leave QC totals untouched", async () => withReport(async page => {
  const quality = await page.locator("#quality").innerHTML();
  await selectFilterField(page, "Condition");
  await valueCheckbox(page, "Mock").uncheck();
  assert.deepEqual(await sampleIds(page), ["001", "002", "020", "030"]);
  await page.locator("#taxonomy-plot-field").selectOption("Latitude");
  assert.deepEqual(await sampleIds(page), ["001", "002", "030", "020"]);
  assert.equal(await valueCheckbox(page, "Mock").isChecked(), false);
  assert.equal(await page.locator("#quality").innerHTML(), quality);
}));

test("multiple metadata filters persist and combine as AND across fields", async () => withReport(async page => {
  await selectFilterField(page, "Condition");
  await page.locator("#taxonomy-filter-none").click();
  await valueCheckbox(page, "Environmental").check();
  assert.deepEqual(await sampleIds(page), ["001", "002", "030"]);
  await selectFilterField(page, "Temperature");
  await page.locator("#taxonomy-filter-none").click();
  await valueCheckbox(page, "20.5").check();
  assert.deepEqual(await sampleIds(page), ["001", "030"]);
  assert.match(await page.locator("#taxonomy-filter-summary").innerText(), /2 active filter/);
  await selectFilterField(page, "Condition");
  assert.equal(await valueCheckbox(page, "Environmental").isChecked(), true);
  assert.equal(await valueCheckbox(page, "Mock").isChecked(), false);
  await page.locator("#taxonomy-rank").selectOption("Division");
  assert.deepEqual(await sampleIds(page), ["001", "030"]);
  await page.locator("#taxonomy-filter-all").click();
  assert.deepEqual(await sampleIds(page), ["001", "030"]); // Temperature remains restricted.
}));

test("SampleID selection preserves leading zeros, empty selections, and reset", async () => withReport(async page => {
  await selectFilterField(page, "SampleID");
  await page.locator("#taxonomy-filter-none").click();
  assert.deepEqual(await sampleIds(page), []);
  assert.match(await page.locator("#taxonomy-explorer-chart").innerText(), /No positive-abundance/);
  await valueCheckbox(page, "001").check();
  assert.deepEqual(await sampleIds(page), ["001"]);
  await valueCheckbox(page, "030").check();
  assert.deepEqual(await sampleIds(page), ["001", "030"]);
  await page.locator("#taxonomy-filter-reset").click();
  assert.deepEqual(await sampleIds(page), ["001", "002", "010", "020", "030"]);
  assert.match(await page.locator("#taxonomy-filter-summary").innerText(), /all samples/);
}));

test("blank metadata is an explicit selectable value, never treated as zero", async () => withReport(async page => {
  await selectFilterField(page, "Latitude");
  await page.locator("#taxonomy-filter-none").click();
  await valueCheckbox(page, "").check();
  assert.deepEqual(await sampleIds(page), ["010", "020"]);
  assert.match(await page.locator("#taxonomy-metadata-values").innerText(), /\(blank\)/);
}));

test("all taxa appear in the legend and table, colours repeat and stay stable", async () => withReport(async page => {
  assert.equal(await page.locator("#taxonomy-explorer-legend .legend-item").count(), 41);
  assert.equal(await page.locator("#taxonomy-explorer-body tr").count(), 41);
  const rare = page.locator("#taxonomy-explorer-legend .legend-item").filter({hasText: /^Taxon 40$/});
  assert.equal(await rare.count(), 1);
  const colours = await page.locator("#taxonomy-explorer-legend .legend-swatch").evaluateAll(swatches => swatches.map(swatch => swatch.style.background));
  assert.equal(new Set(colours).size, 15);
  const before = await rare.locator(".legend-swatch").getAttribute("style");
  await selectFilterField(page, "Condition");
  await valueCheckbox(page, "Mock").uncheck();
  assert.equal(await rare.locator(".legend-swatch").getAttribute("style"), before);
  assert.equal(await page.locator("#taxonomy-explorer-body tr").filter({hasText: /Other/}).count(), 0);
  const denominatorChecks = await page.locator(".taxonomy-stacked-bar").evaluateAll(bars => bars.map(bar => ({
    total: Array.from(bar.children).reduce((sum, segment) => sum + parseFloat(segment.style.height), 0),
    bounded: Array.from(bar.children).every(segment => segment.getBoundingClientRect().top >= bar.getBoundingClientRect().top - 1),
    minHeight: getComputedStyle(bar.firstElementChild).minHeight,
  })));
  for (const bar of denominatorChecks) {
    // Chrome serializes CSS percentages at limited precision (here 100.00001%).
    // Allow less than 0.001 percentage point; payload count totals are checked separately.
    assert.ok(Math.abs(bar.total - 100) < 0.001, `Serialized browser percentages: ${JSON.stringify(bar)}`);
    assert.equal(bar.bounded, true);
    assert.equal(bar.minHeight, "0px");
  }
}));

test("database-specific rank labels keep their original option values", async () => withReport(async page => {
  for (const rank of ["Kingdom", "Phylum"]) {
    assert.equal(await page.locator(`#taxonomy-rank option[value="${rank}"]`).innerText(), `${rank} (SILVA only)`);
  }
  for (const rank of ["Supergroup", "Division", "Subdivision"]) {
    assert.equal(await page.locator(`#taxonomy-rank option[value="${rank}"]`).innerText(), `${rank} (PR2 only)`);
  }
  assert.match(await page.locator('#taxonomy-rank option[value="ProPortal_ASV_Ecotype"]').innerText(), /ProPortal only/);
}));

test("legend scrolls vertically and sample labels remain below the bars", async () => withReport(async page => {
  const layout = await page.evaluate(() => {
    const legend = document.querySelector("#taxonomy-explorer-legend");
    legend.scrollTop = 100;
    return {
      legendHeight: legend.clientHeight,
      scrollHeight: legend.scrollHeight,
      scrollTop: legend.scrollTop,
      labelsBelow: Array.from(document.querySelectorAll(".taxonomy-sample-column")).every(column =>
        column.querySelector(".taxonomy-sample-label").getBoundingClientRect().top >= column.querySelector(".taxonomy-stacked-bar").getBoundingClientRect().bottom),
    };
  });
  assert.ok(layout.legendHeight <= 240);
  assert.ok(layout.scrollHeight > layout.legendHeight);
  assert.ok(layout.scrollTop > 0);
  assert.equal(layout.labelsBelow, true);
  if (process.env.REPORT_SCREENSHOT_DIR) {
    mkdirSync(process.env.REPORT_SCREENSHOT_DIR, {recursive: true});
    // Hide only the navigation in section captures; a sticky bar otherwise overlays tall screenshots.
    await page.locator("nav").evaluate(nav => { nav.style.visibility = "hidden"; });
    await page.locator("#taxa").screenshot({path: join(process.env.REPORT_SCREENSHOT_DIR, "taxonomy-desktop.png")});
    await page.locator("#quality").screenshot({path: join(process.env.REPORT_SCREENSHOT_DIR, "quality-summary.png")});
  }
}));

test("narrow layout keeps its y-axis fixed while bars scroll horizontally", async () => withReport(async page => {
  await page.locator("#taxonomy-plot-field").selectOption("Latitude");
  const checks = await page.evaluate(() => {
    const axis = document.querySelector(".taxonomy-y-axis");
    const scroll = document.querySelector(".taxonomy-scroll");
    const before = axis.getBoundingClientRect().x;
    const title = document.querySelector(".taxonomy-x-axis-title");
    const titleBefore = title.getBoundingClientRect().x;
    scroll.scrollLeft = 120;
    return {
      before, after: axis.getBoundingClientRect().x, moved: scroll.scrollLeft,
      titleBefore, titleAfter: title.getBoundingClientRect().x, title: title.textContent,
      pageWidth: document.documentElement.scrollWidth, viewport: innerWidth,
    };
  });
  assert.equal(checks.before, checks.after);
  assert.equal(checks.titleBefore, checks.titleAfter);
  assert.equal(checks.title, "Latitude");
  assert.ok(checks.moved > 0);
  assert.ok(checks.pageWidth <= checks.viewport);
  if (process.env.REPORT_SCREENSHOT_DIR) {
    // This screenshot-only adjustment does not alter the generated report or its controls.
    await page.locator("nav").evaluate(nav => { nav.style.visibility = "hidden"; });
    await page.locator("#taxa").screenshot({path: join(process.env.REPORT_SCREENSHOT_DIR, "taxonomy-mobile.png")});
  }
}, {viewport: {width: 390, height: 844}}));

test("empty taxonomy shows a clear empty state without browser errors", async () => withReport(async page => {
  const empty = {metadataFields: ["SampleID"], ranks: [], rankLabels: {}, samples: {}};
  const content = readFileSync(fixture, "utf8").replace(
    /(<script type="application\/json" id="taxonomy-explorer-data">)[\s\S]*?(<\/script>)/,
    (_, start, end) => start + JSON.stringify(empty) + end,
  );
  await page.setContent(content);
  assert.deepEqual(await sampleIds(page), []);
  assert.match(await page.locator("#taxonomy-filter-summary").innerText(), /0 of 0 sample/);
  assert.match(await page.locator("#taxonomy-explorer-chart").innerText(), /No positive-abundance/);
}));

test("metadata labels are rendered as text, not executable HTML", async () => withReport(async page => {
  const content = readFileSync(fixture, "utf8").replace(
    /(<script type="application\/json" id="taxonomy-explorer-data">)([\s\S]*?)(<\/script>)/,
    (_, start, encoded, end) => {
      const data = JSON.parse(encoded);
      data.samples["001"].metadata.Temperature = '<img src=x onerror="window.injected=true">';
      return start + JSON.stringify(data).replaceAll("</", "<\\/") + end;
    },
  );
  await page.setContent(content);
  await selectFilterField(page, "Temperature");
  assert.match(await page.locator("#taxonomy-metadata-values").innerText(), /<img src=x/);
  assert.equal(await page.locator("#taxonomy-metadata-values img").count(), 0);
  assert.equal(await page.evaluate(() => window.injected), undefined);
  await page.locator("#taxonomy-plot-field").selectOption("Temperature");
  assert.equal(await page.locator(".taxonomy-sample-label img").count(), 0);
  assert.match((await page.locator(".taxonomy-sample-label").allTextContents()).join(" "), /<img src=x/);
  assert.equal(await page.evaluate(() => window.injected), undefined);
}));

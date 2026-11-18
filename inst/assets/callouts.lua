-- Give the recurring "Exercises" and "Summary" callouts their own look.
-- Each such callout is wrapped in a div (.oma-exercises, .oma-summary)
-- that book.scss styles. Chapters keep using plain callout blocks:
--
--   ::: {.callout-tip}
--   ## Exercises
--   ...
--   :::
--
-- Runs before Quarto parses callouts (at: pre-ast in _quarto.yml), because
-- Quarto drops extra classes from callouts and loses the content of
-- callouts wrapped later in the pipeline.
local kinds = { Exercises = "oma-exercises", Summary = "oma-summary" }

local function is_callout(el)
  for _, c in ipairs(el.classes) do
    if c:match("^callout%-") then return true end
  end
  return false
end

function Div(el)
  if not is_callout(el) then return nil end
  local title = el.attributes["title"]
  if not title then
    local first = el.content[1]
    if first and first.t == "Header" then
      title = pandoc.utils.stringify(first)
    end
  end
  local cls = title and kinds[title]
  if cls then
    return pandoc.Div({ el }, pandoc.Attr("", { cls }))
  end
end

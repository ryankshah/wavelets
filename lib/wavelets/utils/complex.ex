defmodule Wavelets.Utils.Complex do
  @moduledoc """
  Utilities for complex number calculations
  """

  @type t :: {float, float}

  def add({re1, im1}, {re2, im2}), do: {re1 + re2, im1 + im2}
  def sub({re1, im1}, {re2, im2}), do: {re1 - re2, im1 - im2}

  def mul({re1, im1}, {re2, im2}) do
    {re1 * re2 - im1 * im2, re1 * im2 + im1 * re2}
  end

  def scale({re, im}, factor), do: {re * factor, im * factor}
  def conj({re, im}), do: {re, -im}
  def abs({re, im}), do: :math.sqrt(re * re + im * im)
  def phase({re, im}), do: :math.atan2(im, re)
end
